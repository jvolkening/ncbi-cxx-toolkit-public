/* This is a United States Government Work, freely available for public use. */

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <objects/blastdb/Blast_def_line.hpp>
#include <objects/blastdb/Blast_def_line_set.hpp>
#include <objects/blastdb/Blast_db_mask_info.hpp>
#include <objtools/blast/seqdb_reader/seqdb.hpp>
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>
#include <objtools/blast/seqdb_reader/seqdbcommon.hpp>
#include <objtools/blast/seqdb_writer/build_db.hpp>

#include "taxon_exact_cluster.hpp"

USING_NCBI_SCOPE;
USING_SCOPE(objects);

class CRepresentativeSource : public IRawSequenceSource {
public:
    CRepresentativeSource(CSeqDBExpert& db, const vector<int>& oids,
                          CBuildDatabase& writer)
        : m_Db(db), m_Oids(oids)
    {
        m_Db.ListColumns(m_Names);
        vector<string> volumes;
        m_Db.FindVolumePaths(volumes);
        if (volumes.size() > 1 && !m_Names.empty()) {
            NCBI_THROW(CSeqDBException, eArgErr,
                       "Multi-volume databases with columns are not supported");
        }
        for (const auto& name : m_Names) {
            int id = m_Db.GetColumnId(name);
            if (id < 0) {
                NCBI_THROW(CSeqDBException, eFileErr, "Cannot read source column");
            }
            m_ColumnIds.push_back(id);
        }
        m_Db.GetAvailableMaskAlgorithms(m_MaskIds);
        sort(m_MaskIds.begin(), m_MaskIds.end());
        int mask_column = m_Db.GetColumnId("BlastDb/MaskData");
        if (mask_column >= 0 &&
            m_Db.GetColumnMetaData(mask_column).size() != m_MaskIds.size()) {
            NCBI_THROW(CSeqDBException, eArgErr,
                       "Mask algorithm metadata is ambiguous");
        }
        for (int id : m_MaskIds) {
            if (mask_column < 0) {
                NCBI_THROW(CSeqDBException, eArgErr,
                           "Mask metadata is not available as a column");
            }
            const auto& metadata = m_Db.GetColumnMetaData(mask_column);
            auto it = metadata.find(NStr::IntToString(id));
            if (it == metadata.end()) {
                NCBI_THROW(CSeqDBException, eArgErr,
                           "Cannot map source mask algorithm metadata");
            }
            EBlast_filter_program program;
            string name, options;
            m_Db.GetMaskAlgorithmDetails(id, program, name, options);
            int new_id;
            if (program == eBlast_filter_program_other) {
                vector<string> fields;
                NStr::Split(it->second, ":", fields);
                if (fields.size() != 4) {
                    NCBI_THROW(CSeqDBException, eArgErr,
                               "Unsupported custom mask algorithm metadata");
                }
                string description = NStr::Replace(fields[3], string(1, '\1'), ":");
                new_id = writer.RegisterMaskingAlgorithm(name, description, options);
            } else {
                new_id = writer.RegisterMaskingAlgorithm(
                    program, options, name + NStr::IntToString(id));
            }
            m_MaskIdMap[id] = new_id;
            m_MaskMetadata[id] = it->second;
        }
    }

    ~CRepresentativeSource() override
    {
        if (m_Sequence) m_Db.RetSequence(&m_Sequence);
    }

    bool GetNext(CTempString& sequence, CTempString& ambiguities,
                 CRef<CBlast_def_line_set>& deflines,
                 vector<SBlastDbMaskData>& masks, vector<int>& column_ids,
                 vector<CTempString>& column_blobs) override
    {
        if (m_Sequence) {
            m_Db.RetSequence(&m_Sequence);
            m_Sequence = nullptr;
        }
        if (m_Next == m_Oids.size()) return false;
        int oid = m_Oids[m_Next++];
        int found = oid;
        if (!m_Db.CheckOrFindOID(found) || found != oid) {
            NCBI_THROW(CSeqDBException, eFileErr, "Representative OID disappeared");
        }
        deflines = m_Db.GetHdr(oid);
        if (deflines.Empty() || !deflines->CanGet() || deflines->Get().empty()) {
            NCBI_THROW(CSeqDBException, eArgErr, "Representative has no deflines");
        }
        for (const auto& defline : deflines->Get()) {
            if (defline.Empty() || !defline->CanGetSeqid() ||
                defline->GetSeqid().empty()) {
                NCBI_THROW(CSeqDBException, eArgErr,
                           "Representative defline has no sequence ID");
            }
            for (const auto& id : defline->GetSeqid()) {
                if (id.Empty()) {
                    NCBI_THROW(CSeqDBException, eArgErr, "Empty sequence ID");
                }
            }
            if (defline->IsSetOther_info() &&
                defline->GetOther_info().empty()) {
                NCBI_THROW(CSeqDBException, eArgErr,
                           "Representative has an empty other-info field");
            }
        }
        int slength = 0, alength = 0;
        m_Db.GetRawSeqAndAmbig(oid, &m_Sequence, &slength, &alength);
        if (!m_Sequence || slength <= 0 || alength < 0) {
            NCBI_THROW(CSeqDBException, eArgErr, "Representative has no raw sequence");
        }
        sequence = CTempString(m_Sequence, slength);
        ambiguities = CTempString(m_Sequence + slength, alength);
        for (int id : m_MaskIds) {
            CSeqDB::TSequenceRanges ranges;
            m_Db.GetMaskData(oid, id, ranges);
            SBlastDbMaskData data;
            data.algorithm_id = m_MaskIdMap.at(id);
            for (const auto& range : ranges) {
                data.offsets.emplace_back(range.first, range.second);
            }
            masks.push_back(data);
        }
        column_ids = m_ColumnIds;
        m_Blobs.resize(column_ids.size());
        column_blobs.resize(column_ids.size());
        for (size_t i = 0; i < column_ids.size(); ++i) {
            m_Db.GetColumnBlob(column_ids[i], oid, m_Blobs[i]);
            column_blobs[i] = m_Blobs[i].Str();
        }
        return true;
    }

    void GetColumnNames(vector<string>& names) override { names = m_Names; }
    int GetColumnId(const string& name) override { return m_Db.GetColumnId(name); }
    const map<string, string>& GetColumnMetaData(int id) override {
        return m_Db.GetColumnMetaData(id);
    }
    const map<int, int>& GetMaskIdMap() const { return m_MaskIdMap; }
    const map<int, string>& GetMaskMetadata() const { return m_MaskMetadata; }

private:
    CSeqDBExpert& m_Db;
    const vector<int>& m_Oids;
    size_t m_Next = 0;
    const char* m_Sequence = nullptr;
    vector<string> m_Names;
    vector<int> m_ColumnIds, m_MaskIds;
    vector<CBlastDbBlob> m_Blobs;
    map<int, int> m_MaskIdMap;
    map<int, string> m_MaskMetadata;
};

static void s_Verify(CSeqDBExpert& source, const string& output,
                     const vector<int>& representatives,
                     CRepresentativeSource& raw)
{
    CSeqDBExpert result(output, CSeqDB::eNucleotide);
    if (result.GetBlastDbVersion() != eBDB_Version5 ||
        result.GetNumOIDs() != static_cast<int>(representatives.size())) {
        NCBI_THROW(CSeqDBException, eFileErr, "Output database count/version mismatch");
    }
    vector<string> names;
    raw.GetColumnNames(names);
    for (const auto& name : names) {
        if (name == "BlastDb/MaskData") continue;
        int src_id = source.GetColumnId(name);
        int dst_id = result.GetColumnId(name);
        if (dst_id < 0 || source.GetColumnMetaData(src_id) != result.GetColumnMetaData(dst_id)) {
            NCBI_THROW(CSeqDBException, eFileErr, "Output column metadata mismatch: " + name);
        }
    }
    int mask_column = result.GetColumnId("BlastDb/MaskData");
    for (const auto& item : raw.GetMaskIdMap()) {
        if (mask_column < 0) {
            NCBI_THROW(CSeqDBException, eFileErr, "Output mask metadata is missing");
        }
        const auto& metadata = result.GetColumnMetaData(mask_column);
        auto it = metadata.find(NStr::IntToString(item.second));
        if (it == metadata.end() || it->second != raw.GetMaskMetadata().at(item.first)) {
            NCBI_THROW(CSeqDBException, eFileErr, "Output mask metadata mismatch");
        }
    }
    for (size_t i = 0; i < representatives.size(); ++i) {
        int src_oid = representatives[i], dst_oid = static_cast<int>(i), found = dst_oid;
        if (!result.CheckOrFindOID(found) || found != dst_oid) {
            NCBI_THROW(CSeqDBException, eFileErr, "Output OID missing");
        }
        string a, b;
        source.GetSequenceAsString(src_oid, a);
        result.GetSequenceAsString(dst_oid, b);
        CRef<CBlast_def_line_set> ah = source.GetHdr(src_oid), bh = result.GetHdr(dst_oid);
        if (a != b || ah.Empty() || bh.Empty() || !ah->Equals(*bh)) {
            NCBI_THROW(CSeqDBException, eFileErr, "Output sequence or deflines mismatch");
        }
        for (const auto& item : raw.GetMaskIdMap()) {
            CSeqDB::TSequenceRanges ar, br;
            source.GetMaskData(src_oid, item.first, ar);
            result.GetMaskData(dst_oid, item.second, br);
            if (ar.size() != br.size()) {
                NCBI_THROW(CSeqDBException, eFileErr, "Output mask ranges mismatch");
            }
            for (size_t j = 0; j < ar.size(); ++j) {
                if (ar[j].first != br[j].first || ar[j].second != br[j].second) {
                    NCBI_THROW(CSeqDBException, eFileErr, "Output mask ranges mismatch");
                }
            }
        }
        for (const auto& name : names) {
            if (name == "BlastDb/MaskData") continue;
            CBlastDbBlob ab, bb;
            source.GetColumnBlob(source.GetColumnId(name), src_oid, ab);
            result.GetColumnBlob(result.GetColumnId(name), dst_oid, bb);
            if (ab.Str() != bb.Str()) {
                NCBI_THROW(CSeqDBException, eFileErr, "Output column data mismatch: " + name);
            }
        }
    }
}

class CTaxonClusterApp : public CNcbiApplication {
    void Init() override;
    int Run() override;
};

void CTaxonClusterApp::Init()
{
    unique_ptr<CArgDescriptions> args(new CArgDescriptions);
    args->SetUsageContext(GetArguments().GetProgramBasename(),
                          "Cluster nucleotide sequences by taxon and sequence similarity");
    args->AddKey("db", "dbname", "Source BLAST v5 nucleotide database",
                 CArgDescriptions::eString);
    args->AddOptionalKey("out_dir", "directory",
                         "Create a new directory containing a verified v5 database and report",
                         CArgDescriptions::eString);
    args->AddFlag("near",
                  "Experimentally allow up to 1% global edit distance for sequences of 100-10000 bases");
    SetupArgDescriptions(args.release());
}

int CTaxonClusterApp::Run()
{
    CSeqDBExpert db(GetArgs()["db"].AsString(), CSeqDB::eNucleotide);
    if (db.GetBlastDbVersion() != eBDB_Version5) {
        NCBI_THROW(CSeqDBException, eArgErr, "A BLAST v5 database is required");
    }

    vector<int> algorithms;
    db.GetAvailableMaskAlgorithms(algorithms);
    sort(algorithms.begin(), algorithms.end());
    vector<string> volumes;
    db.FindVolumePaths(volumes);
    int mask_column = db.GetColumnId("BlastDb/MaskData");
    if ((!algorithms.empty() && (mask_column < 0 || volumes.size() > 1)) ||
        (mask_column >= 0 &&
         db.GetColumnMetaData(mask_column).size() != algorithms.size())) {
        NCBI_THROW(CSeqDBException, eArgErr,
                   "Cannot unambiguously enumerate source mask algorithms");
    }
    if (mask_column >= 0) {
        const auto& metadata = db.GetColumnMetaData(mask_column);
        for (int id : algorithms) {
            if (metadata.find(NStr::IntToString(id)) == metadata.end()) {
                NCBI_THROW(CSeqDBException, eArgErr,
                           "Cannot map source mask algorithm metadata");
            }
        }
    }

    auto read = [&](int oid) {
        taxon_exact::SRecord record;
        db.GetSequenceAsString(oid, record.sequence);
        CRef<CBlast_def_line_set> headers = db.GetHdr(oid);
        record.complete_taxonomy = headers.NotEmpty() && headers->CanGet() &&
                                   !headers->Get().empty();
        if (record.complete_taxonomy) {
            for (const auto& defline : headers->Get()) {
                if (defline.Empty() || !defline->CanGetSeqid() ||
                    defline->GetSeqid().empty() || !defline->IsSetTaxid() ||
                    TAX_ID_TO(int, defline->GetTaxid()) <= 0) {
                    record.complete_taxonomy = false;
                    break;
                }
                for (const auto& id : defline->GetSeqid()) {
                    if (id.Empty()) record.complete_taxonomy = false;
                }
                record.taxids.push_back(TAX_ID_TO(int, defline->GetTaxid()));
                if (defline->IsSetLinks()) {
                    for (int taxid : defline->GetLinks()) {
                        if (taxid <= 0) record.complete_taxonomy = false;
                        record.taxids.push_back(taxid);
                    }
                }
            }
        }
        for (int algorithm : algorithms) {
            CSeqDB::TSequenceRanges ranges;
            db.GetMaskData(oid, algorithm, ranges);
            taxon_exact::TMasks::value_type mask;
            for (const auto& range : ranges) {
                mask.emplace_back(range.first, range.second);
            }
            record.masks.push_back(std::move(mask));
        }
        return record;
    };

    unique_ptr<CNcbiOfstream> report;
    string out_dir, stage;
    if (GetArgs()["out_dir"]) {
        out_dir = CDirEntry::CreateAbsolutePath(GetArgs()["out_dir"].AsString());
        if (!CDir(out_dir).Create(CDir::fCreate_ErrorIfExists | CDir::fCreate_PermByUmask)) {
            NCBI_THROW(CSeqDBException, eFileErr, "Output directory must not already exist");
        }
        stage = CDirEntry::ConcatPath(out_dir, ".incomplete");
        if (!CDir(stage).Create(CDir::fCreate_ErrorIfExists | CDir::fCreate_PermByUmask)) {
            NCBI_THROW(CSeqDBException, eFileErr, "Cannot create staging directory");
        }
        report.reset(new CNcbiOfstream(CDirEntry::ConcatPath(stage, "clusters.tsv").c_str()));
        if (!*report) NCBI_THROW(CSeqDBException, eFileErr, "Cannot open cluster report");
    }
    CNcbiOstream& out = report ? *report : NcbiCout;
    taxon_exact::CClusterer clusterer(read, std::hash<std::string>(),
                                     GetArgs()["near"].AsBoolean());
    vector<int> representatives;
    out << "representative_oid\tmember_oid\n";
    for (int oid = 0; db.CheckOrFindOID(oid); ++oid) {
        int representative = clusterer.Add(oid, read(oid));
        if (report && representative == oid) representatives.push_back(oid);
        out << representative << '\t' << oid << '\n';
        if (!out) {
            NCBI_THROW(CSeqDBException, eFileErr, "Error writing cluster report");
        }
    }
    out.flush();
    if (!out) {
        NCBI_THROW(CSeqDBException, eFileErr, "Error writing cluster report");
    }
    if (report) {
        report->close();
        if (!*report) NCBI_THROW(CSeqDBException, eFileErr, "Error closing cluster report");
        if (representatives.empty()) {
            NCBI_THROW(CSeqDBException, eArgErr, "Input has no included sequences");
        }
        string prefix = CDirEntry::ConcatPath(stage, "db");
        {
            CBuildDatabase writer(prefix, "Taxon cluster representatives", false,
                                  false, true, false, &NcbiCerr, true,
                                  eBDB_Version5, false, 0);
            writer.SetPreserveRawHeaders();
            CRef<CRepresentativeSource> raw(new CRepresentativeSource(db, representatives, writer));
            if (!writer.AddSequences(*raw) || !writer.EndBuild()) {
                NCBI_THROW(CSeqDBException, eFileErr, "Failed to build output database");
            }
            s_Verify(db, prefix, representatives, *raw);
        }
        if (!CDir(stage).Rename(CDirEntry::ConcatPath(out_dir, "data"))) {
            NCBI_THROW(CSeqDBException, eFileErr, "Cannot publish output directory");
        }
    }
    return 0;
}

int main(int argc, const char* argv[])
{
    return CTaxonClusterApp().AppMain(argc, argv);
}
