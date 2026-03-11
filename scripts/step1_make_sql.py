import collections
import os
import re
import sqlite3
import sys

import pandas as pd
import uuid_utils as uuid
from nanoid import generate

genome = "Human"
assembly = "hg19"

rdf_view_id = str(uuid.uuid7())


datasets = [
    {
        "name": "73 primary",
        "short_name": "73primary",
        "size": 73,
        "institution": "Columbia",
        "index": 1,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "20 ICG",
        "short_name": "20icg",
        "size": 20,
        "institution": "Columbia",
        "index": 2,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "93 Discovery",
        "short_name": "93discovery",
        "size": 93,
        "institution": "Columbia",
        "index": 3,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "29 cell lines",
        "short_name": "29cl",
        "size": 29,
        "institution": "Columbia",
        "index": 4,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "BCCA 150 primary 2024",
        "short_name": "bcca2024-16se",
        "size": 150,
        "institution": "BCCA",
        "index": 5,
        "public_id": str(uuid.uuid7()),
    },
]


HUMAN_CHRS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM",
]

MOUSE_CHRS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chrX",
    "chrY",
    "chrM",
]

CHR_MAP = {
    "Human": {chr: idx + 1 for idx, chr in enumerate(HUMAN_CHRS)},
    "Mouse": {chr: idx + 1 for idx, chr in enumerate(MOUSE_CHRS)},
}

official_symbols = {"human": {}, "mouse": {}}

gene_ids = {"human": [], "mouse": []}
gene_id_map = {"human": {}, "mouse": {}}
prev_gene_id_map = {"human": {}, "mouse": {}}
alias_gene_id_map = {"human": {}, "mouse": {}}

metadata_map = {}

# gene_db_map = {}

file = (
    "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/references/hugo/hugo_20240524.tsv"
)
df_hugo = pd.read_csv(file, sep="\t", header=0, keep_default_na=False)

gene_index = 1

official_symbols = {"human": {}, "mouse": {}}
gene_id_map = {"human": {}, "mouse": {}}

for i, gene_symbol in enumerate(df_hugo["Approved symbol"].values):

    # genes = [gene_id] + list(
    #     filter(
    #         lambda x: x != "",
    #         [x.strip() for x in df_hugo["Previous symbols"].values[i].split(",")],
    #     )
    # )

    hugo = df_hugo["HGNC ID"].values[i]
    ensembl = df_hugo["Ensembl gene ID"].values[i].split(".")[0]
    refseq = df_hugo["RefSeq IDs"].values[i].replace(" ", "")
    ncbi = df_hugo["NCBI Gene ID"].values[i].replace(" ", "")

    info = {
        "index": gene_index,
        "gene_id": hugo,
        "gene_symbol": gene_symbol,
        "ensembl": ensembl,
        "refseq": refseq,
        "ncbi": ncbi,
    }

    official_symbols["human"][hugo] = info
    gene_id_map["human"][gene_symbol.lower()] = hugo
    gene_index += 1


dir = f"../data/modules/wgs"

file = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/wgs/data/human/rdf/hg19/mutation_database/DLBCL_Master_121625.xlsx"

df_samples = pd.read_excel(
    file, sheet_name="Study_Panel_283", header=0, keep_default_na=False
)

# force sample id to be string
df_samples["Sample ID"] = df_samples["Sample ID"].astype(str)

# sample_map = {
#     sample: {"index": i + 1, "public_id": str(uuid.uuid7())}
#     for i, sample in enumerate(df_samples["Sample ID"].values)
# }


# for i, dataset in enumerate(datasets):
#     df_ds = df_samples[df_samples["Dataset"] == dataset]
#     institution = list(df_ds["Institution"].unique())[0]

#     dataset_map[dataset]["index"] = i + 1
#     dataset_map[dataset]["public_id"] = str(uuid.uuid7())


# sampleIdMap = {sample: i for i, sample in enumerate(df_samples["Sample"].values)}

metadata = ["COO", "LymphGen class", "Paired normal DNA", "Institution", "Type"]

metadata_json_map = {
    "COO": {"camelCase": "coo", "short": "c"},
    "LymphGen class": {"camelCase": "lymphgenClass", "short": "lg"},
    "Paired normal DNA": {"camelCase": "pairedNormalDna", "short": "pnd"},
    "Institution": {"camelCase": "institution", "short": "i"},
    "Type": {"camelCase": "type", "short": "t"},
}


metadata_map = {meta: mi + 1 for mi, meta in enumerate(metadata)}


# dataset_dir = os.path.join(dir, shortName)

# if not os.path.exists(dataset_dir):
#    os.makedirs(dataset_dir)

db = os.path.join(dir, f"wgs-20260311.db")

print(db)

if os.path.exists(db):
    os.remove(db)

conn = sqlite3.connect(db)
conn.row_factory = sqlite3.Row
cursor = conn.cursor()


cursor.execute("PRAGMA journal_mode = WAL;")

cursor.execute(
    f"""
    CREATE TABLE genomes (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        name TEXT NOT NULL,
        scientific_name TEXT NOT NULL,
        UNIQUE(name, scientific_name));
    """,
)
cursor.execute("CREATE INDEX idx_genomes_name ON genomes(LOWER(name));")
cursor.execute(
    f"INSERT INTO genomes (id, public_id, name, scientific_name) VALUES (1, '{uuid.uuid7()}', 'Human', 'Homo sapiens');"
)
cursor.execute(
    f"INSERT INTO genomes (id, public_id, name, scientific_name) VALUES (2, '{uuid.uuid7()}', 'Mouse', 'Mus musculus');"
)

genome_map = {
    "human": 1,
    "mouse": 2,
}

cursor.execute(
    f"""
    CREATE TABLE assemblies (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        genome_id INTEGER NOT NULL,
        name TEXT NOT NULL UNIQUE,
        FOREIGN KEY (genome_id) REFERENCES genomes(id) ON DELETE CASCADE);
    """,
)
cursor.execute("CREATE INDEX idx_assemblies_name ON assemblies(LOWER(name));")
cursor.execute("CREATE INDEX idx_assemblies_genome_id ON assemblies (genome_id);")

cursor.execute(
    f"INSERT INTO assemblies (id, public_id, genome_id, name) VALUES (1, '{uuid.uuid7()}', 1, 'hg19');"
)
cursor.execute(
    f"INSERT INTO assemblies (id, public_id, genome_id, name) VALUES (2, '{uuid.uuid7()}', 1, 'GRCh38');"
)
cursor.execute(
    f"INSERT INTO assemblies (id, public_id, genome_id, name) VALUES (3, '{uuid.uuid7()}', 2, 'GRCm39');"
)

assembly_map = {"hg19": 1, "GRCh38": 2, "GRCm39": 3}


cursor.execute(
    f"""
    CREATE TABLE genes (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        genome_id INTEGER NOT NULL,
        gene_id TEXT NOT NULL,
        ensembl TEXT NOT NULL DEFAULT '',
        refseq TEXT NOT NULL DEFAULT '',
        ncbi INTEGER NOT NULL DEFAULT 0,
        gene_symbol TEXT NOT NULL DEFAULT '',
        FOREIGN KEY(genome_id) REFERENCES genomes(id));
    """,
)

genomes = ["human"]


for si, g in enumerate(genomes):
    genome_id = genome_map[g]
    for id in sorted(official_symbols[g.lower()]):
        d = official_symbols[g.lower()][id]

        cursor.execute(
            f"INSERT INTO genes (id, public_id, genome_id, gene_id, ensembl, refseq, ncbi, gene_symbol) VALUES (:id, :public_id, :genome_id, :gene_id, :ensembl, :refseq, :ncbi, :gene_symbol);",
            (
                {
                    "id": d["index"],
                    "public_id": str(uuid.uuid7()),
                    "genome_id": genome_id,
                    "gene_id": d["gene_id"],
                    "ensembl": d["ensembl"],
                    "refseq": d["refseq"],
                    "ncbi": d["ncbi"],
                    "gene_symbol": d["gene_symbol"],
                }
            ),
        )

cursor.execute("COMMIT;")


cursor.execute(
    f"""
    CREATE TABLE institutions (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        name TEXT NOT NULL UNIQUE);
    """,
)
cursor.execute("CREATE INDEX idx_institutions_name ON institutions(LOWER(name));")

cursor.execute(
    f"INSERT INTO institutions (id, public_id, name) VALUES (1, '{uuid.uuid7()}', 'Columbia');"
)
cursor.execute(
    f"INSERT INTO institutions (id, public_id, name) VALUES (2, '{uuid.uuid7()}', 'BCCA');"
)

institution_map = {"columbia": 1, "bcca": 2}


cursor.execute(
    f"""
    CREATE TABLE variant_types (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        name TEXT NOT NULL UNIQUE);
    """,
)
cursor.execute("CREATE INDEX idx_variant_types_name ON variant_types(LOWER(name));")

cursor.execute(
    f"INSERT INTO variant_types (id, public_id, name) VALUES (100, '{uuid.uuid7()}', 'DEL');"
)

cursor.execute(
    f"INSERT INTO variant_types (id, public_id, name) VALUES (200, '{uuid.uuid7()}', 'SNV');"
)

# cursor.execute(
#     f"INSERT INTO variant_types (id, public_id, name) VALUES (300, '{uuid.uuid7()}', 'MNV');"
# )

# it must be last so allow other records to be added before it
cursor.execute(
    f"INSERT INTO variant_types (id, public_id, name) VALUES (300, '{uuid.uuid7()}', 'INS');"
)

# deletions should be ranked before SNVs and insertions for display
# purposes
variant_type_map = {"DEL": 100, "SNV": 200, "INS": 300}


cursor.execute(
    f"""CREATE TABLE info (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL,
    version TEXT NOT NULL DEFAULT ""
    );
    """
)

cursor.execute(
    f"""INSERT INTO info (id, public_id, name, version) VALUES (1, '{uuid.uuid7()}', 'mutations', '1.0.0');"""
)

cursor.execute(
    f"""CREATE TABLE chromosomes (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    genome_id INTEGER NOT NULL,
    chr_id INTEGER NOT NULL,
    name TEXT NOT NULL,
    FOREIGN KEY (genome_id) REFERENCES genomes(id) ON DELETE CASCADE);
"""
)

cursor.execute("CREATE INDEX idx_chromosomes_genome_id ON chromosomes (genome_id);")

chr_map = {
    "Human": {},
    "Mouse": {},
}

chr_index = 1

for chr in HUMAN_CHRS:
    cursor.execute(
        f"INSERT INTO chromosomes (public_id, genome_id, chr_id, name) VALUES ('{str(uuid.uuid7())}', 1, {CHR_MAP['Human'][chr]}, '{chr}');",
    )
    chr_map["Human"][chr] = chr_index
    chr_index += 1

for chr in MOUSE_CHRS:
    cursor.execute(
        f"INSERT INTO chromosomes (public_id, genome_id, chr_id, name) VALUES ('{str(uuid.uuid7())}', 2, {CHR_MAP['Mouse'][chr]}, '{chr}');",
    )
    chr_map["Mouse"][chr] = chr_index
    chr_index += 1

cursor.execute(
    f"""CREATE TABLE permissions (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL
    );
"""
)


cursor.execute(
    f"INSERT INTO permissions (id, public_id, name) VALUES (1, '{rdf_view_id}', 'rdf:view');",
)

cursor.execute(
    f"""CREATE TABLE datasets (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    assembly_id INTEGER NOT NULL,
    institution_id INTEGER NOT NULL,
    name TEXT NOT NULL,
    short_name TEXT NOT NULL,
    samples INTEGER NOT NULL DEFAULT 0,
    mutations INTEGER NOT NULL DEFAULT 0,
    description TEXT NOT NULL DEFAULT "",
    FOREIGN KEY(assembly_id) REFERENCES assemblies(id),
    FOREIGN KEY(institution_id) REFERENCES institutions(id)
    );
"""
)

cursor.execute("CREATE INDEX idx_datasets_name ON datasets (LOWER(name));")
cursor.execute("CREATE INDEX idx_datasets_assembly_id ON datasets (assembly_id);")
cursor.execute("CREATE INDEX idx_datasets_institution_id ON datasets (institution_id);")

cursor.execute(
    f"""CREATE TABLE dataset_permissions (
    dataset_id INTEGER NOT NULL,
    permission_id INTEGER NOT NULL,
    PRIMARY KEY (dataset_id, permission_id),
    FOREIGN KEY(dataset_id) REFERENCES datasets(id),
    FOREIGN KEY(permission_id) REFERENCES permissions(id)
    );
"""
)

cursor.execute(
    "CREATE INDEX idx_dataset_permissions_dataset_id ON dataset_permissions (dataset_id);"
)
cursor.execute(
    "CREATE INDEX idx_dataset_permissions_permission_id ON dataset_permissions (permission_id);"
)


cursor.execute(
    f""" CREATE TABLE samples (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL,
    coo TEXT NOT NULL DEFAULT "",
    lymphgen_class TEXT NOT NULL DEFAULT "",
    paired_normal_dna TEXT NOT NULL DEFAULT "",
    type TEXT NOT NULL DEFAULT ""
    );
    """
)

# not required as unique constraint is on public_id creates index
# cursor.execute("CREATE INDEX idx_samples_public_id ON samples (public_id);")
cursor.execute("CREATE INDEX idx_samples_name ON samples (LOWER(name));")

cursor.execute(
    f""" CREATE TABLE sample_metadata (
    sample_id INTEGER NOT NULL,
    metadata_id INTEGER NOT NULL,
    value TEXT NOT NULL,
    PRIMARY KEY (sample_id, metadata_id),
    FOREIGN KEY(sample_id) REFERENCES samples(id),
    FOREIGN KEY(metadata_id) REFERENCES metadata(id)
    );
    """
)

cursor.execute(
    "CREATE INDEX idx_sample_metadata_sample_id ON sample_metadata (sample_id);"
)
cursor.execute(
    "CREATE INDEX idx_sample_metadata_metadata_id ON sample_metadata (metadata_id);"
)


# a sample can be in multiple datasets, so we need a many-to-many relationship
cursor.execute(
    f""" CREATE TABLE dataset_samples (
    dataset_id INTEGER NOT NULL,
    sample_id INTEGER NOT NULL,
    PRIMARY KEY (dataset_id, sample_id),
    FOREIGN KEY(dataset_id) REFERENCES datasets(id),
    FOREIGN KEY(sample_id) REFERENCES samples(id)
    );
    """
)

cursor.execute(
    "CREATE INDEX idx_dataset_samples_dataset_id ON dataset_samples (dataset_id);"
)
cursor.execute(
    "CREATE INDEX idx_dataset_samples_sample_id ON dataset_samples (sample_id);"
)


cursor.execute(
    f""" CREATE TABLE variants (
    id INTEGER PRIMARY KEY,
    chr_id INTEGER NOT NULL,
    gene_id INTEGER,
    variant_type_id INTEGER NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    ref TEXT NOT NULL,
    tum TEXT NOT NULL,
    hgvs_c               TEXT NOT NULL DEFAULT "",               -- coding HGVS notation
    hgvs_p               TEXT NOT NULL DEFAULT "",               -- protein HGVS notation
    consequence          TEXT NOT NULL DEFAULT "",               -- missense, frameshift, etc.
    clinical_significance TEXT NOT NULL DEFAULT "",              -- Pathogenic, VUS, Benign...
    FOREIGN KEY(chr_id) REFERENCES chromosomes(id),
    FOREIGN KEY(gene_id) REFERENCES genes(id),
    FOREIGN KEY(variant_type_id) REFERENCES variant_types(id)
    );
    """
)

cursor.execute(
    f"""CREATE INDEX idx_variants_chr_id_start_end ON variants (chr_id, start, end);"""
)

cursor.execute(
    f"""CREATE INDEX idx_variants_variant_type_id ON variants (variant_type_id);"""
)

cursor.execute(f"""CREATE INDEX idx_variants_gene_id ON variants (gene_id);""")


cursor.execute(
    f"""CREATE TABLE annotation (
  id             INTEGER  PRIMARY KEY,
  variant_id     INTEGER  NOT NULL REFERENCES variants(id),
  source         TEXT     NOT NULL,  -- ClinVar, gnomAD, SIFT...
  source_version TEXT,               -- 2024-01, v4.0...
  key            TEXT     NOT NULL,  -- e.g. CLNSIG, sift_score
  value          TEXT     NOT NULL,

  UNIQUE (variant_id, source, key)
);"""
)

cursor.execute(
    f"""CREATE INDEX idx_annotation_variant ON annotation (variant_id, LOWER(source));"""
)

cursor.execute(
    f"""CREATE INDEX idx_annotation_variant_source_key ON annotation (variant_id, LOWER(source), LOWER(key));"""
)

cursor.execute(
    f""" CREATE TABLE sample_variants (
    dataset_id INTEGER NOT NULL,
    sample_id INTEGER NOT NULL,
    variant_id INTEGER NOT NULL,
    t_alt_count INTEGER NOT NULL DEFAULT -1,
    t_depth INTEGER NOT NULL DEFAULT -1,
    vaf FLOAT NOT NULL DEFAULT -1,
    PRIMARY KEY (dataset_id, sample_id, variant_id),
    FOREIGN KEY(dataset_id) REFERENCES datasets(id),
    FOREIGN KEY(sample_id) REFERENCES samples(id),
    FOREIGN KEY(variant_id) REFERENCES variants(id)
    );
    """
)


# cursor.execute(
#     f"""CREATE INDEX idx_mutations_gene_symbol ON mutations (LOWER(hugo_gene_symbol)); """
# )

# cursor.execute(
#     f""" CREATE TABLE sample_mutations (
#     sample_id INTEGER NOT NULL,
#     mutation_id INTEGER NOT NULL,
#     PRIMARY KEY (sample_id, mutation_id),
#     FOREIGN KEY(sample_id) REFERENCES samples(id),
#     FOREIGN KEY(mutation_id) REFERENCES mutations(id)
#     );
#     """
# )


# for meta in metadata:
#     uuid = uuid7()
#     cursor.execute(
#         f"INSERT INTO metadata (id, public_id, name, short_name) VALUES ({metadata_map[meta]}, '{public_id}', '{meta}', '{metadata_json_map[meta]["camelCase"]}');",
#     )

genome = "Human"
# file = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/wgs/data/human/rdf/hg19/mutation_database/bcca2024_73primary_29cl_20icg_hg19/bcca2024_73primary_29cl_20icg_hg19.maf.txt"
file = "/home/antony/development/ngs/wgs/bcca2024/bcca2024-16se_73primary_29cl_20icg_hg19.fixed.maf.txt"

df = pd.read_csv(file, sep="\t", header=0, keep_default_na=False)
df["Sample"] = df["Sample"].astype(str)

# sort by Chromosome col
# df = df.sort_values(by="Chromosome")

variant_map = {}
sample_map = {}

for di, dataset in enumerate(datasets):
    dataset_index = dataset["index"]
    dataset_id = dataset["public_id"]
    institution = dataset["institution"]
    name = dataset["name"]
    short_name = dataset["short_name"]
    size = dataset["size"]

    if institution.lower() not in institution_map:
        idx = len(institution_map) + 1
        institution_map[institution.lower()] = idx
        cursor.execute(
            f"INSERT INTO institutions (id, public_id, name) VALUES ({idx}, '{str(uuid.uuid7())}', '{institution}');",
        )

    institution_id = institution_map[institution.lower()]

    dfd = df[df["Dataset"] == short_name]

    mutation_counts = dfd.shape[0]

    cursor.execute(
        f"INSERT INTO datasets (id, public_id, assembly_id, institution_id, name, short_name, samples, mutations) VALUES ({dataset_index}, '{dataset_id}', {assembly_map[assembly]}, {institution_id}, '{name}', '{short_name}', {size}, {mutation_counts});",
    )

    cursor.execute(
        f"INSERT INTO dataset_permissions (dataset_id, permission_id) VALUES ({dataset_index}, 1);",
    )

    for i, row in dfd.iterrows():
        # mutation_uuid = str(uuid.uuid7())
        # generate("0123456789abcdefghijklmnopqrstuvwxyz", 12)

        sample = row["Sample"]
        # remove _allele_1 etc from sample name
        sample = re.sub(r"_.+", "", sample)

        if sample not in sample_map:
            sample_map[sample] = len(sample_map) + 1

            sample_index = sample_map[sample]
            sample_id = str(uuid.uuid7())

            df_samples_d = df_samples[df_samples["Sample ID"] == sample]
            print(sample, df_samples_d)

            coo = df_samples_d["COO class"].values[0]

            if "nd" in coo.lower():
                coo = "NA"

            lymphgen = df_samples_d["LymphGen class"].values[0]
            paired = df_samples_d["Paired normal DNA"].values[0]
            # ins = df_samples["Institution"].values[i]
            sample_type = df_samples_d["Sample type"].values[0]

            cursor.execute(
                f"INSERT INTO samples (id, public_id, name, coo, lymphgen_class, paired_normal_dna, type) VALUES ({sample_index}, '{sample_id}', '{sample}', '{coo}', '{lymphgen}', '{paired}', '{sample_type}') ON CONFLICT DO NOTHING;",
            )

        sample_index = sample_map[sample]

        # map samples to datasets

        cursor.execute(
            f"INSERT INTO dataset_samples (dataset_id, sample_id) VALUES ({dataset_index}, {sample_index}) ON CONFLICT DO NOTHING;",
        )

        if dataset_index == 1 or dataset_index == 2:
            # sample is also added to the 93 discovery dataset if part of the 73 or 20 icg
            cursor.execute(
                f"INSERT INTO dataset_samples (dataset_id, sample_id) VALUES (3, {sample_index}) ON CONFLICT DO NOTHING;",
            )

        chr = row["Chromosome"]

        chr = chr.upper().replace("MT", "M")
        chr = chr.replace("CHR", "")
        chr = "chr" + chr

        if chr not in chr_map[genome]:
            continue

        chr_id = chr_map[genome][chr]

        start = row["Start_Position"]
        end = row["End_Position"]
        ref = row["Reference_Allele"]
        tum = row["Tumor_Seq_Allele2"]
        gene = row["Hugo_Symbol"]
        vaf = row["VAF"]
        db = row["Dataset"]
        hgvs_c = row["DNAChange"]
        hgvs_p = row["AAChange"]
        consequence = row["Variant_Classification"]

        # default to assuming no gene found
        gene_id = "NULL"

        for g in gene.split(";"):
            if g != "":
                id = gene_id_map[genome.lower()].get(g.lower(), -1)

                if id != -1:
                    # print(f"Gene: {g}, Gene ID: {gene_id}", "dataset", dataset)

                    # change gene id from null to the actual gene id, if found
                    gene_id = official_symbols[genome.lower()][id]["index"]
                    break

        if vaf == "na":
            vaf = -1

        # variant_type = dfd["Variant_Type"].values[i]

        t_alt_count = row["t_alt_count"]
        t_depth = row["t_depth"]

        # if t_alt_count == "na":
        #    t_alt_count = -1

        # if t_depth == "na":
        #    t_depth = -1

        snvs = []

        # split concatenated snps
        if len(ref) > 1 and len(tum) > 1 and len(ref) == len(tum):
            for idx in range(len(ref)):
                ref_i = ref[idx]
                tum_i = tum[idx]

                if ref_i == tum_i:
                    continue

                start_i = start + idx
                end_i = end - (len(ref) - idx - 1)

                snvs.append(
                    {"start": start_i, "end": end_i, "ref": ref_i, "tum": tum_i}
                )

            variant_type = "SNV"
        else:
            snvs.append({"start": start, "end": end, "ref": ref, "tum": tum})

            if ref[0] == "-":
                variant_type = "INS"
            elif tum[0] == "-":
                variant_type = "DEL"
            else:
                # print(f"Unknown variant type for ref: {ref}, tum: {tum} {start}")
                variant_type = "SNV"
                # sys.exit(1)

        if variant_type not in variant_type_map:
            print(variant_type)
            variant_type_id = len(variant_type_map) + 1
            variant_type_map[variant_type] = variant_type_id
            cursor.execute(
                f"INSERT INTO variant_types (id, public_id, name) VALUES ({variant_type_id}, '{uuid.uuid7()}', '{variant_type}');"
            )

        variant_type_id = variant_type_map[variant_type]

        for snv in snvs:
            # we unique mutations by their genomic location and change to reduce repeats
            variant_key = "|".join(
                [
                    str(chr_id),
                    str(variant_type_id),
                    str(gene_id),
                    str(snv["start"]),
                    str(snv["end"]),
                    snv["ref"],
                    snv["tum"],
                    hgvs_p,
                    consequence,
                ]
            )

            if variant_key not in variant_map:
                mutation_index = len(variant_map) + 1
                variant_map[variant_key] = mutation_index

                cursor.execute(
                    f"""INSERT INTO variants (id, chr_id, variant_type_id, gene_id, start, end, ref, tum, hgvs_c, hgvs_p, consequence) VALUES (
                    {mutation_index}, 
                    {chr_id}, 
                    {variant_type_id}, 
                    {gene_id}, 
                    {snv['start']}, 
                    {snv['end']}, 
                    '{snv['ref']}', 
                    '{snv['tum']}', 
                    '{hgvs_c}', 
                    '{hgvs_p}', 
                    '{consequence}');
                    """
                )

            mutation_index = variant_map[variant_key]

            # so we can merge mutations from different tables, use the public_id as foreign key

            # if we are in the 73 primary or 20 ICG datasets, also add to the 93 discovery dataset
            if dataset_index == 1 or dataset_index == 2:
                datasets = [dataset_index, 3]
            else:
                # only add variant to the current dataset
                datasets = [dataset_index]

            for dsi in datasets:
                cursor.execute(
                    f"""INSERT INTO sample_variants (dataset_id, sample_id, variant_id, t_alt_count, t_depth, vaf) VALUES (
                    {dsi}, 
                    {sample_index}, 
                    {mutation_index}, 
                    {t_alt_count}, 
                    {t_depth}, 
                    {vaf}) ON CONFLICT DO NOTHING;"""
                )


conn.commit()
conn.close()
