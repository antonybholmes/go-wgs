import collections
import os
import sqlite3
import sys

import pandas as pd
import uuid_utils as uuid
from nanoid import generate

genome = "Human"
assembly = "hg19"

rdf_view_id = str(uuid.uuid7())


idMap = {"20_icg": "20icg", "29_cell_lines": "29cl", "73_bcca": "73primary"}

renameMap = {
    "20_icg": "20 ICG",
    "29_cell_lines": "29 cell lines",
    "73_bcca": "73 primary",
}

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

dir = f"../data/modules/mutations"

file = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/wgs/data/human/rdf/hg19/mutation_database/93primary_29cl_dlbcl_hg19/samples.txt"

df_samples = pd.read_csv(file, sep="\t", header=0, keep_default_na=False)
datasets = list(sorted(df_samples["Dataset"].unique()))
sample_map = {
    sample: {"index": i + 1, "public_id": str(uuid.uuid7())}
    for i, sample in enumerate(df_samples["Sample"].values)
}

dataset_map = {}

for i, dataset in enumerate(datasets):
    df_ds = df_samples[df_samples["Dataset"] == dataset]
    institution = list(df_ds["Institution"].unique())[0]

    dataset_map[dataset] = {
        "index": i + 1,
        "public_id": str(uuid.uuid7()),
        "institution": institution,
    }

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

db = os.path.join(dir, f"mutations.db")

print(db)

if os.path.exists(db):
    os.remove(db)

conn = sqlite3.connect(db)
conn.row_factory = sqlite3.Row
cursor = conn.cursor()

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
    f"""CREATE TABLE metadata (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL,
    version TEXT NOT NULL DEFAULT ""
    );
    """
)

cursor.execute(
    f"""INSERT INTO metadata (id, public_id, name, version) VALUES (1, '{uuid.uuid7()}', 'mutations', '1.0.0');"""
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
    institution TEXT NOT NULL DEFAULT "",
    name TEXT NOT NULL,
    short_name TEXT NOT NULL,
    mutations INTEGER NOT NULL DEFAULT 0,
    description TEXT NOT NULL DEFAULT "",
    FOREIGN KEY(assembly_id) REFERENCES assemblies(id)
    );
"""
)

cursor.execute("CREATE INDEX idx_datasets_name ON datasets (LOWER(name));")
cursor.execute("CREATE INDEX idx_datasets_assembly_id ON datasets (assembly_id);")

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
    dataset_id INTEGER NOT NULL,
    name TEXT NOT NULL,
    coo TEXT NOT NULL DEFAULT "",
    lymphgen_class TEXT NOT NULL DEFAULT "",
    paired_normal_dna TEXT NOT NULL DEFAULT "",
    type TEXT NOT NULL DEFAULT "",
    FOREIGN KEY(dataset_id) REFERENCES datasets(id)
    );
    """
)
cursor.execute("CREATE INDEX idx_samples_dataset_id ON samples (dataset_id);")
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

cursor.execute(
    f""" CREATE TABLE mutations (
    id INTEGER PRIMARY KEY,
    sample_id INTEGER NOT NULL,
    hugo_gene_symbol TEXT NOT NULL DEFAULT '',
    variant_classification TEXT NOT NULL DEFAULT '',
    variant_type TEXT NOT NULL DEFAULT '',
    chr_id INTEGER NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    ref TEXT NOT NULL,
    tum TEXT NOT NULL,
    t_alt_count INTEGER NOT NULL DEFAULT -1,
    t_depth INTEGER NOT NULL DEFAULT -1,
    vaf FLOAT NOT NULL DEFAULT -1,
    FOREIGN KEY(sample_id) REFERENCES samples(id),
    FOREIGN KEY(chr_id) REFERENCES chromosomes(id)
    );
    """
)
cursor.execute(
    f"""CREATE INDEX idx_mutations_chr_id_start_end ON mutations (chr_id, start, end);"""
)
cursor.execute(
    f"""CREATE INDEX idx_mutations_gene_symbol ON mutations (LOWER(hugo_gene_symbol)); """
)


# for meta in metadata:
#     uuid = uuid7()
#     cursor.execute(
#         f"INSERT INTO metadata (id, public_id, name, short_name) VALUES ({metadata_map[meta]}, '{public_id}', '{meta}', '{metadata_json_map[meta]["camelCase"]}');",
#     )

genome = "Human"
file = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/wgs/data/human/rdf/hg19/mutation_database/93primary_29cl_dlbcl_hg19/93primary_29cl_rename_samples_hg19.maf.txt"

df = pd.read_csv(file, sep="\t", header=0, keep_default_na=False)

# sort by Chromosome col
# df = df.sort_values(by="Chromosome")

for di, dataset in enumerate(datasets):
    dataset_index = dataset_map[dataset]["index"]
    dataset_id = dataset_map[dataset]["public_id"]
    institution = dataset_map[dataset]["institution"]

    short_name = idMap[dataset]
    df_samples_d = df_samples[df_samples["Dataset"] == dataset]

    dfd = df[df["Sample"].isin(df_samples_d["Sample"])]

    print(short_name)

    mutation_counts = dfd.shape[0]

    name = renameMap[dataset]

    cursor.execute(
        f"INSERT INTO datasets (id, public_id, assembly_id, institution, name, short_name, mutations) VALUES ({dataset_index}, '{dataset_id}', {assembly_map[assembly]}, '{institution}', '{name}', '{short_name}', {mutation_counts});",
    )

    cursor.execute(
        f"INSERT INTO dataset_permissions (dataset_id, permission_id) VALUES ({dataset_index}, 1);",
    )

    for i in range(df_samples_d.shape[0]):
        sample = df_samples_d["Sample"].values[i]
        sample_index = sample_map[sample]["index"]
        sample_id = sample_map[sample]["public_id"]

        coo = df_samples_d["COO"].values[i]

        if "nd" in coo.lower():
            coo = "NA"

        lymphgen = df_samples_d["LymphGen class"].values[i]
        paired = df_samples_d["Paired normal DNA"].values[i]
        # ins = df_samples["Institution"].values[i]
        sample_type = df_samples_d["Sample type"].values[i]

        cursor.execute(
            f"INSERT INTO samples (id, public_id, dataset_id, name, coo, lymphgen_class, paired_normal_dna, type) VALUES ({sample_index}, '{sample_id}', {dataset_index}, '{sample}', '{coo}', '{lymphgen}', '{paired}', '{sample_type}') ON CONFLICT DO NOTHING;",
        )

    for i in range(dfd.shape[0]):
        mutation_uuid = str(
            uuid.uuid7()
        )  # generate("0123456789abcdefghijklmnopqrstuvwxyz", 12)

        chr = dfd["Chromosome"].values[i]

        chr = chr.upper().replace("MT", "M")
        chr = chr.replace("CHR", "")
        chr = "chr" + chr

        if chr not in chr_map[genome]:
            continue

        chr_id = chr_map[genome][chr]

        start = dfd["Start_Position"].values[i]
        end = dfd["End_Position"].values[i]
        ref = dfd["Reference_Allele"].values[i]
        tum = dfd["Tumor_Seq_Allele2"].values[i]
        vaf = dfd["VAF"].values[i]
        db = dfd["Database"].values[i]

        if vaf == "na":
            vaf = -1

        variant_type = dfd["Variant_Type"].values[i]

        t_alt_count = dfd["t_alt_count"].values[i]
        t_depth = dfd["t_depth"].values[i]

        if t_alt_count == "na":
            t_alt_count = -1

        if t_depth == "na":
            t_depth = -1

        # so we can merge mutations from different tables, use the public_id as foreign key
        cursor.execute(
            f"INSERT INTO mutations (sample_id, chr_id, start, end, ref, tum, t_alt_count, t_depth, variant_type, vaf) VALUES ({sample_index}, {chr_id}, {start}, {end}, '{ref}', '{tum}', {t_alt_count}, {t_depth}, '{variant_type}', {vaf});",
        )


conn.commit()
conn.close()
