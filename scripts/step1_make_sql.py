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
        "description": "Variants in 73 primary cases.",
        # "size": 73,
        "institution": "Columbia",
        "index": 1,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "20 ICG",
        "short_name": "20icg",
        "description": "Variants in 20 ICG cases.",
        # "size": 20,
        "institution": "Columbia",
        "index": 2,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "93 Discovery",
        "short_name": "93discovery",
        "description": "The 73primary + 20icg cases used in the Bal Nature paper.",
        # "size": 93,
        "institution": "Columbia",
        "index": 3,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "29 cell lines",
        "short_name": "29cl-dlbcl",
        "description": "Variants in 29 DLBCL cell lines.",
        # "size": 29,
        "institution": "Columbia",
        "index": 4,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "BCCA 2024 16 SE",
        "short_name": "bcca2024-16se",
        "description": "Variant data from 16 SE regions in the BCCA 2024 dataset.",
        # "size": 222,
        "institution": "BCCA",
        "index": 5,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "WGS 122",
        "short_name": "122wgs",
        "description": "The 73primary + 20icg + 29cl cases.",
        # "size": 122,
        "institution": "Columbia",
        "index": 6,
        "public_id": str(uuid.uuid7()),
    },
    {
        "name": "BCCA 2024 16 SE Unique",
        "short_name": "bcca2024-16se-unique",
        "description": "Primary samples from the BCCA 2024 dataset in 16 SE regions that are not in the 93 discovery dataset.",
        # "size": 150,
        "institution": "BCCA",
        "index": 7,
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

file = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/wgs/data/human/rdf/hg19/mutation_database/DLBCL_Master_040726.xlsx"

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

db = os.path.join(dir, f"wgs-20260415.db")

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


# cursor.execute(
#     f"""
#     CREATE TABLE features (
#         id INTEGER PRIMARY KEY,
#         public_id TEXT NOT NULL UNIQUE,
#         name TEXT NOT NULL UNIQUE)
#     """,
# )
# cursor.execute("CREATE INDEX idx_features_name ON features(LOWER(name));")

# cursor.execute(
#     f"INSERT INTO features (id, public_id, name) VALUES (1, '{uuid.uuid7()}', 'gene');"
# )
# cursor.execute(
#     f"INSERT INTO features (id, public_id, name) VALUES (2, '{uuid.uuid7()}', 'transcript');"
# )
# cursor.execute(
#     f"INSERT INTO features (id, public_id, name) VALUES (3, '{uuid.uuid7()}', 'exon');"
# )
# cursor.execute(
#     f"INSERT INTO features (id, public_id, name) VALUES (4, '{uuid.uuid7()}', 'CDS');"
# )

# feature_map = {"gene": 1, "transcript": 2, "exon": 3, "CDS": 4}

cursor.execute(
    f"""
    CREATE TABLE biotypes (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        name TEXT NOT NULL UNIQUE);
    """
)

cursor.execute("CREATE INDEX idx_biotypes_name ON biotypes(LOWER(name));")

# assume this is a hugo gene symbol. We can add no hugo
# symbols later
cursor.execute(
    f"""
    CREATE TABLE genes (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        genome_id INTEGER NOT NULL,
        biotype_id INTEGER NOT NULL,
        gene_id TEXT NOT NULL,
        gene_symbol TEXT NOT NULL DEFAULT '',
        is_hugo BOOLEAN NOT NULL DEFAULT 0,
        is_canonical BOOLEAN NOT NULL DEFAULT 0,
        FOREIGN KEY(genome_id) REFERENCES genomes(id),
        FOREIGN KEY(biotype_id) REFERENCES biotypes(id));
    """
)

cursor.execute("CREATE INDEX idx_genes_gene_id ON genes(LOWER(gene_id));")
cursor.execute("CREATE INDEX idx_genes_gene_symbol ON genes(LOWER(gene_symbol));")
cursor.execute("CREATE INDEX idx_genes_genome_id ON genes(genome_id);")
cursor.execute("CREATE INDEX idx_genes_biotype_id ON genes(biotype_id);")

cursor.execute(
    f"""
    CREATE TABLE transcripts (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        gene_id INTEGER NOT NULL,
        name TEXT NOT NULL,
        mane_refseq TEXT NOT NULL DEFAULT "",
        mane_status TEXT NOT NULL DEFAULT "",
        exons INTEGER NOT NULL DEFAULT 0,
        FOREIGN KEY(gene_id) REFERENCES genes(id));
    """
)

cursor.execute("CREATE INDEX idx_transcripts_name ON transcripts(LOWER(name));")


df_genes = pd.read_csv(
    "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/references/gencode/grch37/gencode_v48lift37_basic_transcripts.tsv",
    sep="\t",
    header=0,
)

biotype_map = {}
gene_map = {}
transcript_map = {}

for i, row in df_genes.iterrows():
    gene_id = row["gene_id"]
    gene_symbol = row["gene_symbol"]
    hgnc_id = row["hgnc_id"]
    is_hugo = bool(hgnc_id)
    is_canonical = bool(row["is_canonical"])
    transcript = row["transcript_id"]
    exons = row["exons"]
    biotype = row["biotype"]

    if biotype not in biotype_map:
        biotype_map[biotype] = len(biotype_map) + 1
        biotype_db_id = biotype_map[biotype]
        cursor.execute(
            f"INSERT INTO biotypes (id, public_id, name) VALUES (:id, :public_id, :name);",
            (
                {
                    "id": biotype_db_id,
                    "public_id": str(uuid.uuid7()),
                    "name": biotype,
                }
            ),
        )

    biotype_db_id = biotype_map[biotype]

    if gene_id not in gene_map:
        gene_map[gene_id] = len(gene_map) + 1
        gene_db_id = gene_map[gene_id]
        cursor.execute(
            f"INSERT INTO genes (id, public_id, genome_id, biotype_id, gene_id, gene_symbol, is_hugo, is_canonical) VALUES (:id, :public_id, :genome_id, :biotype_id, :gene_id, :gene_symbol, :is_hugo, :is_canonical);",
            (
                {
                    "id": gene_db_id,
                    "public_id": str(uuid.uuid7()),
                    "genome_id": genome_map[genome.lower()],
                    "biotype_id": biotype_db_id,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "is_hugo": is_hugo,
                    "is_canonical": is_canonical,
                }
            ),
        )

    gene_db_id = gene_map[gene_id]

    if transcript not in transcript_map:
        transcript_map[transcript] = len(transcript_map) + 1
        transcript_db_id = transcript_map[transcript]
        cursor.execute(
            f"INSERT INTO transcripts (id, public_id, gene_id, name, exons) VALUES (:id, :public_id, :gene_id, :name, :exons);",
            (
                {
                    "id": transcript_db_id,
                    "public_id": str(uuid.uuid7()),
                    "gene_id": gene_db_id,
                    "name": transcript,
                    "exons": exons,
                }
            ),
        )


genomes = ["human"]


# for si, g in enumerate(genomes):
#     genome_id = genome_map[g]
#     for id in sorted(official_symbols[g.lower()]):
#         d = official_symbols[g.lower()][id]

#         cursor.execute(
#             f"INSERT INTO genes (id, public_id, genome_id, gene_id, gene_symbol, is_hugo, is_canonical) VALUES (:id, :public_id, :genome_id, :gene_id, :gene_symbol, :is_hugo, :is_canonical);",
#             (
#                 {
#                     "id": d["index"],
#                     "public_id": str(uuid.uuid7()),
#                     "genome_id": genome_id,
#                     "gene_id": d["gene_id"],
#                     "gene_symbol": d["gene_symbol"],
#                     "is_hugo": d["is_hugo"],
#                     "is_canonical": d["is_canonical"],
#                 }
#             ),
#         )

# cursor.execute("COMMIT;")


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
    f"""INSERT INTO info (id, public_id, name, version) VALUES (1, '{uuid.uuid7()}', 'wgs', '1.0.0');"""
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

#     clinical_significance TEXT NOT NULL DEFAULT "",              -- Pathogenic, VUS, Benign...
cursor.execute(
    f""" CREATE TABLE variants (
    id INTEGER PRIMARY KEY,
    chr_id INTEGER NOT NULL,
    gene_id INTEGER NOT NULL,
    transcript_id INTEGER NOT NULL,
    variant_type_id INTEGER NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    exon INTEGER NOT NULL DEFAULT 0,
    ref TEXT NOT NULL,
    tum TEXT NOT NULL,
    hgvs_c               TEXT NOT NULL DEFAULT "",               -- coding HGVS notation
    hgvs_p               TEXT NOT NULL DEFAULT "",               -- protein HGVS notation
    consequence          TEXT NOT NULL DEFAULT "",               -- missense, frameshift, etc.
    FOREIGN KEY(chr_id) REFERENCES chromosomes(id),
    FOREIGN KEY(gene_id) REFERENCES genes(id),
    FOREIGN KEY(transcript_id) REFERENCES transcripts(id),
    FOREIGN KEY(variant_type_id) REFERENCES variant_types(id)
    );
    """
)

cursor.execute(
    f""" CREATE TABLE secondary_variants (
    id INTEGER PRIMARY KEY,
    variant_id INTEGER NOT NULL,
    gene_id INTEGER,
    hgvs_c               TEXT NOT NULL DEFAULT "",               -- coding HGVS notation
    hgvs_p               TEXT NOT NULL DEFAULT "",               -- protein HGVS notation
    consequence          TEXT NOT NULL DEFAULT "",               -- missense, frameshift, etc.
    FOREIGN KEY(variant_id) REFERENCES variants(id),
    FOREIGN KEY(gene_id) REFERENCES genes(id)
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
    f"""CREATE INDEX idx_variants_transcript_id ON variants (transcript_id);"""
)


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

# match a given sample to a variant along with observed counts for that sample
cursor.execute(
    f""" CREATE TABLE sample_variants (
    id INTEGER PRIMARY KEY,
    sample_id INTEGER NOT NULL,
    variant_id INTEGER NOT NULL,
    n_alt_count INTEGER NOT NULL DEFAULT 0,
    n_depth INTEGER NOT NULL DEFAULT 0,
    t_alt_count INTEGER NOT NULL DEFAULT 0,
    t_depth INTEGER NOT NULL DEFAULT 0,
    vaf FLOAT NOT NULL DEFAULT 0,
    UNIQUE (sample_id, variant_id),
    FOREIGN KEY(sample_id) REFERENCES samples(id),
    FOREIGN KEY(variant_id) REFERENCES variants(id)
    );
    """
)

# datasets can share variants, so we need a many-to-many relationship
cursor.execute(
    f""" CREATE TABLE dataset_sample_variants (
    dataset_id INTEGER NOT NULL,
    sample_variant_id INTEGER NOT NULL,
    PRIMARY KEY (dataset_id, sample_variant_id),
    FOREIGN KEY(dataset_id) REFERENCES datasets(id),
    FOREIGN KEY(sample_variant_id) REFERENCES sample_variants(id)
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
file = "/home/antony/development/ngs/wgs/bcca2024/v4/bcca2024-16se_73primary_29cl_20icg_hg19.maf.vep.splice.reordered.fix-datasets.v5.txt"  # "/home/antony/development/ngs/wgs/bcca2024/bcca2024-16se_73primary_29cl_20icg_hg19.fixed.maf.txt"


# sort by Chromosome col
# df = df.sort_values(by="Chromosome")


variant_map = {}
sample_map = {}

sample_variant_map = {}


for di, dataset in enumerate(datasets):
    dataset_index = dataset["index"]
    dataset_id = dataset["public_id"]
    institution = dataset["institution"]
    name = dataset["name"]
    short_name = dataset["short_name"]
    # size = dataset["size"]
    description = dataset["description"]

    print(
        f"Processing dataset {dataset_index}: {name} with institution {institution}..."
    )

    if institution.lower() not in institution_map:
        idx = len(institution_map) + 1
        institution_map[institution.lower()] = idx
        cursor.execute(
            f"INSERT INTO institutions (id, public_id, name) VALUES ({idx}, '{str(uuid.uuid7())}', '{institution}');",
        )

    institution_id = institution_map[institution.lower()]

    cursor.execute(
        f"""INSERT INTO datasets (id, public_id, assembly_id, institution_id, name, short_name, description) VALUES (
                {dataset_index}, 
                '{dataset_id}', 
                {assembly_map[assembly]}, 
                {institution_id}, 
                '{name}', 
                '{short_name}', 
                '{description}');
            """,
    )

    cursor.execute(
        f"INSERT INTO dataset_permissions (dataset_id, permission_id) VALUES ({dataset_index}, 1);",
    )

for df in pd.read_csv(
    file, sep="\t", header=0, keep_default_na=False, chunksize=200000
):
    df["Tumor_Sample_Barcode"] = df["Tumor_Sample_Barcode"].astype(str)

    for di, dataset in enumerate(datasets):
        dataset_index = dataset["index"]
        dataset_id = dataset["public_id"]
        institution = dataset["institution"]
        name = dataset["name"]
        short_name = dataset["short_name"]
        # size = dataset["size"]
        description = dataset["description"]

        print(
            f"Processing dataset {dataset_index}: {name} with institution {institution}..."
        )

        institution_id = institution_map[institution.lower()]

        # dfd = df[df["Dataset"].str.split("|").apply(lambda x: short_name in x)]
        dfd = df[
            df["Dataset"].str.contains(rf"(?:^|\|){short_name}(?:\||$)", regex=True)
        ]

        if dfd.shape[0] == 0:
            print(f"No rows found for dataset {name}, skipping...")
            continue

        for i, row in dfd.iterrows():
            # mutation_uuid = str(uuid.uuid7())
            # generate("0123456789abcdefghijklmnopqrstuvwxyz", 12)

            sample = row["Tumor_Sample_Barcode"]
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
                    f"""INSERT INTO samples (id, public_id, name, coo, lymphgen_class, paired_normal_dna, type) VALUES (
                    {sample_index}, 
                    '{sample_id}', 
                    '{sample}', 
                    '{coo}', 
                    '{lymphgen}', 
                    '{paired}', 
                    '{sample_type}') 
                    ON CONFLICT DO NOTHING;
                """,
                )

            sample_index = sample_map[sample]

            # map samples to datasets

            cursor.execute(
                f"INSERT INTO dataset_samples (dataset_id, sample_id) VALUES ({dataset_index}, {sample_index}) ON CONFLICT DO NOTHING;",
            )

            # if dataset_index == 1 or dataset_index == 2:
            #     # sample is also added to the 93 discovery dataset if part of the 73 or 20 icg
            #     cursor.execute(
            #         f"INSERT INTO dataset_samples (dataset_id, sample_id) VALUES (3, {sample_index}) ON CONFLICT DO NOTHING;",
            #     )

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
            gene_id = row["VEP_Gene_ID"]  # row["Entrez_Gene_Id"]
            gene_symbol = row["VEP_Gene_Symbol"]  # row["Hugo_Symbol"]
            is_hugo = row["VEP_Is_Hugo_Gene"]
            is_canonical = row["VEP_Is_Canonical"]
            exon = row["VEP_Exon"]
            exons = row["VEP_Total_Exons"]
            biotype = row["VEP_Biotype"]

            transcript = row["VEP_Transcript"]
            mane_refseq = row["MANE_RefSeq"]
            mane_status = row["MANE_status"]
            vaf = row["VAF"]
            db = row["Dataset"]
            hgvs_c = row["VEP_HGVSc"]  # row["DNAChange"]
            hgvs_p = row["VEP_HGVSp"]  # row["AAChange"]
            consequence = row["VEP_Variant_Classification"]

            # reduce space by using empty string instead of NA or . for missing values
            if consequence == "NA" or consequence == ".":
                consequence = ""

            if hgvs_c == "NA" or hgvs_c == ".":
                hgvs_c = ""

            if hgvs_p == "NA" or hgvs_p == ".":
                hgvs_p = ""

            if gene_symbol == "NA" or gene_symbol == ".":
                gene_symbol = ""

            if transcript == "NA" or transcript == ".":
                transcript = ""

            if mane_refseq == "NA" or mane_refseq == ".":
                mane_refseq = ""

            if mane_status == "NA" or mane_status == ".":
                mane_status = ""

            if biotype == "NA" or biotype == ".":
                biotype = ""

            if exon == "NA" or exon == ".":
                exon = 0

            if exons == "NA" or exons == ".":
                exons = 0

            if biotype not in biotype_map:
                biotype_map[biotype] = len(biotype_map) + 1
                biotype_db_id = biotype_map[biotype]
                cursor.execute(
                    f"INSERT INTO biotypes (id, public_id, name) VALUES (:id, :public_id, :name);",
                    (
                        {
                            "id": biotype_db_id,
                            "public_id": str(uuid.uuid7()),
                            "name": biotype,
                        }
                    ),
                )

            biotype_db_id = biotype_map[biotype]

            if gene_id not in gene_map:
                gene_map[gene_id] = len(gene_map) + 1
                gene_db_id = gene_map[gene_id]
                cursor.execute(
                    f"INSERT INTO genes (id, public_id, genome_id, biotype_id, gene_id, gene_symbol, is_hugo, is_canonical) VALUES (:id, :public_id, :genome_id, :biotype_id, :gene_id, :gene_symbol, :is_hugo, :is_canonical);",
                    (
                        {
                            "id": gene_db_id,
                            "public_id": str(uuid.uuid7()),
                            "genome_id": genome_map[genome.lower()],
                            "biotype_id": biotype_db_id,
                            "gene_id": gene_id,
                            "gene_symbol": gene_symbol,
                            "is_hugo": is_hugo,
                            "is_canonical": is_canonical,
                        }
                    ),
                )

            gene_db_id = gene_map[gene_id]

            if transcript not in transcript_map:
                transcript_map[transcript] = len(transcript_map) + 1
                transcript_db_id = transcript_map[transcript]
                cursor.execute(
                    f"INSERT INTO transcripts (id, public_id, gene_id, name, mane_refseq, mane_status, exons) VALUES (:id, :public_id, :gene_id, :name, :mane_refseq, :mane_status, :exons);",
                    (
                        {
                            "id": transcript_db_id,
                            "public_id": str(uuid.uuid7()),
                            "gene_id": gene_db_id,
                            "name": transcript,
                            "mane_refseq": mane_refseq,
                            "mane_status": mane_status,
                            "exons": exons,
                        }
                    ),
                )

            transcript_db_id = transcript_map[transcript]

            # for g in gene_symbol.split(";"):
            #     if g != "":
            #         id = gene_id_map[genome.lower()].get(g.lower(), -1)

            #         if id != -1:
            #             # print(f"Gene: {g}, Gene ID: {gene_id}", "dataset", dataset)

            #             # change gene id from null to the actual gene id, if found
            #             gene_id = official_symbols[genome.lower()][id]["index"]
            #             break

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
                # print(
                #     gene_id,
                #     transcript,
                #     snv["start"],
                #     snv["end"],
                #     snv["ref"],
                #     snv["tum"],
                #     hgvs_p,
                #     consequence,
                # )

                variant_key = "|".join(
                    [
                        str(chr_id),
                        str(variant_type_id),
                        gene_id,
                        transcript,
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
                        f"""INSERT INTO variants (id, chr_id, variant_type_id, gene_id, transcript_id, start, end, exon, ref, tum, hgvs_c, hgvs_p, consequence) VALUES (
                        {mutation_index}, 
                        {chr_id}, 
                        {variant_type_id}, 
                        {gene_db_id}, 
                        {transcript_db_id}, 
                        {snv['start']}, 
                        {snv['end']}, 
                        {exon},
                        '{snv['ref']}', 
                        '{snv['tum']}', 
                        '{hgvs_c}', 
                        '{hgvs_p}', 
                        '{consequence}');
                        """
                    )

                variant_index = variant_map[variant_key]

                # so we can merge mutations from different tables, use the public_id as foreign key

                # if we are in the 73 primary or 20 ICG datasets, also add to the 93 discovery dataset
                # if dataset_index == 1 or dataset_index == 2:
                #     datasets = [dataset_index, 3]
                # else:
                #     # only add variant to the current dataset
                #     datasets = [dataset_index]

                # for dsi in datasets:

                sample_variant_key = (
                    f"{sample_index}|{variant_index}|{t_alt_count}|{t_depth}"
                )

                if sample_variant_key not in sample_variant_map:
                    sample_variant_map[sample_variant_key] = len(sample_variant_map) + 1
                    sample_variant_index = sample_variant_map[sample_variant_key]
                    cursor.execute(
                        f"""INSERT INTO sample_variants (id, sample_id, variant_id, t_alt_count, t_depth, vaf) VALUES (
                        {sample_variant_index},
                        {sample_index}, 
                        {variant_index}, 
                        {t_alt_count}, 
                        {t_depth}, 
                        {vaf}) ON CONFLICT DO NOTHING;
                        """
                    )

                sample_variant_index = sample_variant_map[sample_variant_key]

                cursor.execute(
                    f"""INSERT INTO dataset_sample_variants (dataset_id, sample_variant_id) VALUES (
                    {dataset_index}, 
                    {sample_variant_index}) ON CONFLICT DO NOTHING;
                    """
                )


conn.commit()
conn.close()
