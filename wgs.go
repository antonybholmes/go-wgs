package wgs

import (
	"database/sql"
	"fmt"

	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-web"
	"github.com/antonybholmes/go-web/auth/sqlite"
	"github.com/antonybholmes/go-web/auth/token"
)

type (
	ExpressionDataIndex struct {
		ProbeIds    []string
		EntrezIds   []string
		GeneSymbols []string
	}

	ExpressionData struct {
		Exp    [][]float64
		Header []string
		Index  ExpressionDataIndex
	}

	MutationsReq struct {
		Assembly string        `json:"assembly"`
		Location *dna.Location `json:"location"`
		Samples  []string      `json:"samples"`
	}

	SampleMetadata struct {
		Name  string `json:"name"`
		Value string `json:"value"`
	}

	Sample struct {
		PublicId        string `json:"id"`
		Name            string `json:"name"`
		Dataset         string `json:"dataset"`
		COO             string `json:"coo,omitempty"`
		LymphgenClass   string `json:"lymphgenClass,omitempty"`
		PairedNormalDna string `json:"pairedNormalDna,omitempty"`
		Type            string `json:"type,omitempty"`
	}

	DatasetStats struct {
		Id       int `json:"id"`
		Samples  int `json:"samples"`
		Variants int `json:"variants"`
	}

	Dataset struct {
		PublicId    string    `json:"id"`
		Genome      string    `json:"genome"`
		Assembly    string    `json:"assembly"`
		Institution string    `json:"institution,omitempty"`
		Name        string    `json:"name"`
		ShortName   string    `json:"shortName"`
		Description string    `json:"description"`
		File        string    `json:"-"`
		Samples     []*Sample `json:"samples"`
		Id          int       `json:"-"`
		Variants    int       `json:"variants"`
	}

	Variant struct {
		Sample      string   `json:"sample"`
		Chr         string   `json:"chr"`
		Ref         string   `json:"ref"`
		Tum         string   `json:"tum"`
		Gene        string   `json:"gene"`
		Transcript  string   `json:"transcript"`
		Type        string   `json:"type"`
		Consequence string   `json:"consequence"`
		HGVSp       string   `json:"hgvsP,omitempty"`
		HGVSc       string   `json:"hgvsC,omitempty"`
		Datasets    []string `json:"datasets"`
		Start       int      `json:"start"`
		End         int      `json:"end"`
		Exon        int      `json:"exon"`
		Exons       int      `json:"exons"`
		NRefCount   int      `json:"nRefCount"`
		NAltCount   int      `json:"nAltCount"`
		NDepth      int      `json:"nDepth"`
		TRefCount   int      `json:"tRefCount"`
		TAltCount   int      `json:"tAltCount"`
		TDepth      int      `json:"tDepth"`
		Id          int      `json:"-"`
		Vaf         float64  `json:"vaf"`
		Y           int      `json:"y,omitempty"`
	}

	DatasetResults struct {
		Dataset  string     `json:"dataset"`
		Variants []*Variant `json:"variants"`
	}

	PileupResults struct {
		Location *dna.Location     `json:"location"`
		Datasets []string          `json:"datasets"`
		Pileup   []*PileupLocation `json:"pileup"`
	}

	PileupLocation struct {
		Variants []*Variant `json:"variants"`
		Start    int        `json:"start"`
		Index    int        `json:"index"`
		// keep track of max Y at position
		Y int `json:"-"`
	}

	DatasetSearchResults struct {
		Dataset  string     `json:"dataset"`
		Variants []*Variant `json:"variants"`
	}

	VariantSearchResults struct {
		Location *dna.Location `json:"location"`
		Datasets []string      `json:"datasets"`
		Variants []*Variant    `json:"variants"`
	}

	MAFLocation struct {
		Start int `json:"start"`
		Index int `json:"index"`
		// keep track of max Y at position
		Count int `json:"count"`
	}

	MAFResults struct {
		Location *dna.Location  `json:"location"`
		Datasets []string       `json:"datasets"`
		MAFs     []*MAFLocation `json:"mafs"`
		Samples  int            `json:"samples"`
		Alleles  int            `json:"alleles"`
	}

	WGSDB struct {
		db   *sql.DB
		file string
	}
)

const (
	MaxVariants = 1000

	DatasetsSql = `SELECT DISTINCT
		d.id,
		d.public_id AS public_id,
		g.name AS genome,
		a.name AS assembly,
		ins.name AS institution,
		d.name,
		d.short_name,
		d.description,
		s.public_id AS sample_id,
		s.name AS sample_name,
		s.coo,
		s.lymphgen_class,
		s.paired_normal_dna
		FROM datasets d
		JOIN assemblies a ON d.assembly_id = a.id
		JOIN genomes g ON a.genome_id = g.id
		JOIN institutions ins ON d.institution_id = ins.id
		JOIN dataset_samples ds ON d.id = ds.dataset_id
		JOIN samples s ON s.id = ds.sample_id
		JOIN dataset_permissions dp ON d.id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE 
			<<PERMISSIONS>>
			AND LOWER(a.name) = :assembly
		ORDER BY
			d.name, 
			s.name`

	DatasetSamplesAndVariantsSql = `SELECT
		d.id,
		COUNT(DISTINCT sv.sample_id) AS samples,
		COUNT(DISTINCT sv.variant_id) AS variants
		FROM datasets d
		JOIN dataset_sample_variants dsv ON d.id = dsv.dataset_id
		JOIN sample_variants sv ON dsv.sample_variant_id = sv.id
		GROUP BY d.id`

	BaseSearchSamplesSql = `SELECT	
		s.id,
		g.name AS genome,
		a.name AS assembly,
		d.name as dataset_name,
		s.name as sample_name,
		FROM samples s
		JOIN datasets d ON s.dataset_id = d.id
		JOIN assemblies a ON d.assembly_id = a.id
		JOIN genomes g ON a.genome_id = g.id
		JOIN dataset_permissions dp ON s.dataset_id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE 
			<<PERMISSIONS>>
			AND LOWER(a.name) = :assembly`

	AllSamplesSql = BaseSearchSamplesSql +
		` ORDER BY
			d.name, 
			s.name`

	SearchSamplesSql = BaseSearchSamplesSql +
		` AND (s.id = :id OR d.id = :id OR d.name LIKE :q OR s.name LIKE :q)
		ORDER BY
			d.name, 
			s.name`

	// SampleMetadataSql = `SELECT
	// 	md.name,
	// 	sm.value
	// 	FROM sample_metadata sm
	// 	JOIN metadata md ON md.id = sm.metadata_id
	// 	WHERE sm.sample_id = :sample_id
	// 	ORDER BY md.name`

	// order by vt.id so that dels are first, then snvs, then ins so that they will be plotted in that order
	FindVariantsSql = `
		SELECT DISTINCT
			v.id as variant_id,
			c.name AS chr,
			v.start, 
			v.end,
			v.ref, 
			v.tum,
			v.hgvs_c,
			v.hgvs_p,
			v.consequence,
			v.exon,
			vt.name AS variant_type,
			g.gene_symbol AS gene_symbol,
			t.name AS transcript,
			t.exons,
			d.public_id AS dataset,
			s.public_id AS sample,
			sv.n_depth - sv.n_alt_count AS n_ref_count, 
			sv.n_alt_count, 
			sv.n_depth,
			sv.t_depth - sv.t_alt_count AS t_ref_count, 
			sv.t_alt_count, 
			sv.t_depth,
			ROUND(sv.vaf, 4) AS vaf
		FROM variants v
		JOIN chromosomes c ON c.id = v.chr_id
		JOIN variant_types vt ON vt.id = v.variant_type_id
		JOIN genes g ON v.gene_id = g.id
		JOIN transcripts t ON v.transcript_id = t.id
		JOIN sample_variants sv ON sv.variant_id = v.id
		JOIN samples s ON sv.sample_id = s.id
		JOIN dataset_sample_variants dsv ON dsv.sample_variant_id = sv.id
		JOIN datasets d ON dsv.dataset_id = d.id
		JOIN dataset_permissions dp ON d.id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE
			<<PERMISSIONS>>
			AND <<DATASETS>>
			AND c.name = :chr AND v.end >= :start AND v.start <= :end
		ORDER BY vt.id, v.start, v.id, s.id, d.id
		LIMIT :limit
	`

	// simpler for counting mutations at each position for maf output
	// since we don't need all the details just the position and count
	MAFSql = `
		SELECT DISTINCT
			sv.sample_id AS sample_id,
			sv.variant_id AS variant_id,
			v.start, 
			v.end
		FROM variants v
		JOIN chromosomes c ON c.id = v.chr_id
		JOIN sample_variants sv ON sv.variant_id = v.id
		JOIN samples s ON sv.sample_id = s.id
		JOIN dataset_sample_variants dsv ON dsv.sample_variant_id = sv.id
		JOIN datasets d ON dsv.dataset_id = d.id
		JOIN dataset_permissions dp ON d.id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE
			<<PERMISSIONS>>
			AND <<DATASETS>>
			AND c.name = :chr AND v.end >= :start AND v.start <= :end
		ORDER BY v.start
		LIMIT :limit
	`

	SampleCountSql = `
		SELECT COUNT(DISTINCT s.id)
		FROM samples s
		JOIN dataset_samples ds ON s.id = ds.sample_id
		JOIN datasets d ON ds.dataset_id = d.id
		JOIN dataset_permissions dp ON d.id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE
			<<PERMISSIONS>>
			AND <<DATASETS>>
	`

	// FindSampleDatasetsSql = `
	// 	SELECT DISTINCT
	// 		d.name
	// 	FROM samples s
	// 	JOIN dataset_samples ds ON s.id = ds.sample_id
	// 	JOIN datasets d ON ds.dataset_id = d.id
	// 	WHERE
	// 		s.id = :sample
	//`
)

func MakeInDatasetsSql(query string, datasetIds []string, namedArgs *[]any) string {
	inPlaceholders := make([]string, len(datasetIds))

	for i, datasetId := range datasetIds {
		ph := fmt.Sprintf("did%d", i+1)
		inPlaceholders[i] = ":" + ph
		*namedArgs = append(*namedArgs, sql.Named(ph, datasetId))
	}

	return strings.Replace(query, "<<DATASETS>>", "d.public_id IN ("+strings.Join(inPlaceholders, ",")+")", 1)
}

func (variant *Variant) Clone() *Variant {
	var ret Variant = Variant{
		Id:        variant.Id,
		Chr:       variant.Chr,
		Start:     variant.Start,
		End:       variant.End,
		Ref:       variant.Ref,
		Tum:       variant.Tum,
		TRefCount: variant.TRefCount,
		TAltCount: variant.TAltCount,
		Gene:      variant.Gene,
		TDepth:    variant.TDepth,
		Type:      variant.Type,
		Vaf:       variant.Vaf,
		Sample:    variant.Sample,
		Datasets:  variant.Datasets,
		Y:         variant.Y,
	}

	return &ret
}

func NewWGSDB(file string) *WGSDB {

	db := sys.Must(sql.Open(sys.Sqlite3DB, file+sys.SqliteDSN))

	return &WGSDB{file: file, db: db}
}

// func (mdb *WGSDB) Dir() string {
// 	return mdb.dir
// }

func (wdb *WGSDB) Datasets(assembly string, isAdmin bool, permissions []string) ([]*Dataset, error) {
	// determine size of all datasets

	rows, err := wdb.db.Query(DatasetSamplesAndVariantsSql)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	// here we get stats for all datasets for speed since we don't
	// want to do separate queries and we also don't want to do
	// the more complex query with permissions since we just want the stats
	mapDatasetInfo := make(map[int]DatasetStats)

	var datasetId int
	var samples int
	var variants int

	for rows.Next() {
		err := rows.Scan(&datasetId, &samples, &variants)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		mapDatasetInfo[datasetId] = DatasetStats{Id: datasetId, Samples: samples, Variants: variants}
	}

	// here only datasets with permissions are collated
	namedArgs := []any{sql.Named("assembly", web.FormatParam(assembly))}

	query := sqlite.MakePermissionsSql(DatasetsSql, isAdmin, permissions, &namedArgs)

	rows, err = wdb.db.Query(query, namedArgs...)

	if err != nil {
		return nil, err //fmt.Errorf("there was an error with the database query")
	}

	defer rows.Close()

	ret := make([]*Dataset, 0, 10)

	var currentDataset *Dataset = nil
	var currentSample *Sample = nil

	for rows.Next() {
		var dataset Dataset
		var sample Sample

		err := rows.Scan(&dataset.Id,
			&dataset.PublicId,
			&dataset.Genome,
			&dataset.Assembly,
			&dataset.Institution,
			&dataset.Name,
			&dataset.ShortName,
			&dataset.Description,
			&sample.PublicId,
			&sample.Name,
			&sample.COO,
			&sample.LymphgenClass,
			&sample.PairedNormalDna,
		)

		if err != nil {
			return nil, err //fmt.Errorf("there was an error with the database records")
		}

		// add the variant count if it exists
		if info, ok := mapDatasetInfo[dataset.Id]; ok {
			dataset.Variants = info.Variants
		}

		if currentDataset == nil || currentDataset.PublicId != dataset.PublicId {
			// if a new dataset, add to the list
			currentDataset = &dataset
			currentDataset.Samples = make([]*Sample, 0, 100)
			ret = append(ret, currentDataset)
		}

		if currentSample == nil || currentSample.PublicId != sample.PublicId {
			sample.Dataset = dataset.PublicId

			// new sample, init metadata map
			//sample.Metadata = make(map[string]string)
			currentSample = &sample
			currentDataset.Samples = append(currentDataset.Samples, currentSample)
			//currentDataset.Samples = append(currentDataset.Samples, &sample)
		}
	}

	return ret, nil
}

func (wdb *WGSDB) Search(assembly string,
	location *dna.Location,
	datasetIds []string,
	isAdmin bool,
	permissions []string) ([]*Variant, error) {

	namedArgs := []any{
		sql.Named("chr", location.Chr()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("limit", MaxVariants)}

	query := sqlite.MakePermissionsSql(FindVariantsSql, isAdmin, permissions, &namedArgs)

	query = MakeInDatasetsSql(query, datasetIds, &namedArgs)

	rows, err := wdb.db.Query(query, namedArgs...)

	//log.Debug().Msgf("querying mutations db with %v %s", isAdmin, query)

	if err != nil {
		//log.Debug().Msgf("error querying mutations db: %s", err)
		return nil, err
	}

	defer rows.Close()

	results := make([]*Variant, 0, 100)

	//var currentDatasetResults *DatasetResults

	//v.id as variant_id,
	// c.name AS chr,
	// v.start,
	// v.end,
	// v.ref,
	// v.tum,
	// v.hgvs_c,
	// v.hgvs_p,
	// v.consequence,
	// vt.name AS variant_type,
	// COALESCE(g.gene_symbol, '') AS gene_symbol,
	// d.id AS dataset_id,
	// d.name AS dataset,
	// s.id as sample_id,
	// s.name AS sample,
	// sv.t_depth - sv.t_alt_count AS ref_count,
	// sv.t_alt_count,
	// sv.t_depth,
	// sv.vaf

	var currentVariant *Variant = nil

	var dataset string

	for rows.Next() {
		var variant Variant

		err := rows.Scan(
			&variant.Id,
			&variant.Chr,
			&variant.Start,
			&variant.End,
			&variant.Ref,
			&variant.Tum,
			&variant.HGVSc,
			&variant.HGVSp,
			&variant.Consequence,
			&variant.Exon,
			&variant.Type,
			&variant.Gene,
			&variant.Transcript,
			&variant.Exons,
			&dataset,
			&variant.Sample,
			&variant.NRefCount,
			&variant.NAltCount,
			&variant.NDepth,
			&variant.TRefCount,
			&variant.TAltCount,
			&variant.TDepth,
			&variant.Vaf,
		)

		//log.Debug().Msgf("found mutation %v", mutation)

		if err != nil {
			return nil, err
		}

		// if there are multiple datasets for the same variant, they will be grouped together
		// by keeping the first variant and adding datasets to it as they come in since ordered by dataset
		// if the variant changes then it must be a new entry however if the variant id stays the same but
		// sample changes then it is also a new entry since same variant can be in multiple samples but
		// we want to keep them separate for plotting purposes.
		if currentVariant == nil || variant.Id != currentVariant.Id || variant.Sample != currentVariant.Sample {
			currentVariant = &variant
			currentVariant.Datasets = make([]string, 0, 5)
			results = append(results, currentVariant)
		}

		// if same variant, just add dataset to list
		currentVariant.Datasets = append(currentVariant.Datasets, dataset)

		// since ordered by dataset we can group here
		// if currentDatasetResults == nil || currentDatasetResults.Dataset != variant.Dataset {
		// 	// create a new stores for this dataset and set to current
		// 	currentDatasetResults = &DatasetResults{Dataset: variant.Dataset, Variants: make([]*Variant, 0, 100)}
		// 	results.DatasetResults = append(results.DatasetResults, currentDatasetResults)
		// }

	}

	return results, nil
}

// MAF returns the minor allele frequency at a location for a set of datasets by counting the number of samples with a mutation at
// each position and dividing by the total number of alleles in those samples (2 per sample for diploid).
// This returns the counts and alleles assuming consumer will translate to MAF as needed.
func (wdb *WGSDB) MAF(assembly string,
	location *dna.Location,
	datasetIds []string,
	isAdmin bool,
	user *token.AuthUserJwtClaims) (*MAFResults, error) {

	namedArgs := []any{}

	query := sqlite.MakePermissionsSql(SampleCountSql, isAdmin, user.Permissions, &namedArgs)

	query = MakeInDatasetsSql(query, datasetIds, &namedArgs)

	var samples int
	err := wdb.db.QueryRow(query, namedArgs...).Scan(&samples)

	namedArgs = []any{
		sql.Named("chr", location.Chr()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End()),
		sql.Named("limit", MaxVariants)}

	query = sqlite.MakePermissionsSql(MAFSql, isAdmin, user.Permissions, &namedArgs)

	query = MakeInDatasetsSql(query, datasetIds, &namedArgs)

	rows, err := wdb.db.Query(query, namedArgs...)

	if err != nil {
		return nil, err
	}

	defer rows.Close()

	var variantId int
	var sampleId int
	var start int
	var end int

	mafMap := make(map[int]*MAFLocation)
	mafs := make([]*MAFLocation, 0, 200)

	for rows.Next() {

		err := rows.Scan(
			&sampleId,
			&variantId,
			&start,
			&end,
		)

		//log.Debug().Msgf("found mutation %v", mutation)

		if err != nil {
			return nil, err
		}

		// if there are multiple datasets for the same variant, they will be grouped together
		// by keeping the first variant and adding datasets to it as they come in since ordered by dataset
		// if the variant changes then it must be a new entry however if the variant id stays the same but
		// sample changes then it is also a new entry since same variant can be in multiple samples but
		// we want to keep them separate for plotting purposes.

		for pos := start; pos <= end; pos++ {
			maf, ok := mafMap[pos]

			// for each position in the variant, add to maf map if not there and then increment count for maf at that position
			if !ok {
				maf = &MAFLocation{Start: pos, Index: pos - location.Start() + 1, Count: 0}
				mafMap[pos] = maf
				mafs = append(mafs, maf)
			}

			maf.Count += 1
		}

	}

	results := MAFResults{Location: location, Datasets: datasetIds, Samples: samples, Alleles: 2 * samples, MAFs: mafs}

	return &results, nil
}

// func (wdb *WGSDB) getDatasetsForSample(sample string) ([]string, error) {

// 	namedArgs := []any{
// 		sql.Named("sample", sample)}

// 	rows, err := wdb.db.Query(FindSampleDatasetsSql, namedArgs...)

// 	if err != nil {
// 		return nil, err
// 	}

// 	defer rows.Close()

// 	results := make([]string, 0, 5)

// 	var dataset string

// 	for rows.Next() {

// 		err := rows.Scan(&dataset)

// 		if err != nil {
// 			return nil, err
// 		}

// 		results = append(results, dataset)
// 	}

// 	return results, nil
// }

// Converts a list of variants at a location into pileup format for plotting
func GetPileup(location *dna.Location, datasets []string, variants []*Variant) (*PileupResults, error) {
	// first lets fix deletions and insertions
	// for _, datasetResults := range search.DatasetResults {
	// 	for _, mutation := range datasetResults.Variants {
	// 		// change for sorting purposes so that ins always comes last
	// 		switch mutation.Type {
	// 		case "INS":
	// 			mutation.Type = "2:INS"
	// 			// modify the output so that is begins with a caret to indicate
	// 			// an insertion
	// 			mutation.Tum = fmt.Sprintf("^%s", mutation.Tum)
	// 		case "DEL":
	// 			mutation.Type = "3:DEL"
	// 		default:
	// 			mutation.Type = "1:SNV"
	// 		}
	// 	}
	// }

	start := location.Start()
	end := location.End()

	// put together by position, type, tum

	//yMap := make(map[int]int)
	//pileupMap := make(map[int]map[string]map[string][]*Variant)

	// init pileups
	pileups := make([]*PileupLocation, location.Len())

	for i := range location.Len() {
		pileups[i] = &PileupLocation{Start: location.Start() + i, Index: i, Variants: make([]*Variant, 0, 100)}
	}

	y := -1

	//the relative position from target start
	pos := -1
	var pileup *PileupLocation

	for _, variant := range variants {
		switch variant.Type {
		case "DEL":
			y = -1
			for i := range variant.End - variant.Start + 1 {
				p := variant.Start + i
				pos = p - start

				if p < start || p > end {
					continue
				}

				pileup = pileups[pos]

				// since we are in position order,
				// find the next height slot for this
				// entry
				if i == 0 {
					pileups[pos], y = nextYSlot(pileup)
				} else {
					// update every position of deletion to have
					// new y of start so that deletion will remain
					// on its own row with no break
					pileups[pos].Y = pileups[pos-1].Y
				}

				pileups[pos] = addToPileup(variant, y, pileup)
			}

		case "INS":
			pos = variant.Start - start
			pileup = pileups[pos]
			pileups[pos], y = nextYSlot(pileup)
			pileups[pos] = addToPileup(variant, y, pileup)
		default:
			// deal with concatenated snps e.g. ACG>TGT
			//tum := []rune(mutation.Tum)
			y = -1
			for i, c := range variant.Tum {
				// clone and change tumor
				mut2 := variant.Clone()
				mut2.Tum = string(c)
				mut2.Type = "SNV"
				p := mut2.Start + i
				pos = p - start

				if p < start || p > end {
					continue
				}

				pileup = pileups[pos]

				if i == 0 {
					pileups[pos], y = nextYSlot(pileup)
				} else {
					// update every position of deletion to have
					// new y of start so that deletion will remain
					// on its own row with no break
					pileups[pos].Y = pileups[pos-1].Y
				}

				pileups[pos] = addToPileup(mut2, y, pileup)
			}

		}

	}

	// keep only positions with something in them
	filtered := make([]*PileupLocation, 0, len(pileups))

	for _, loc := range pileups {
		if len(loc.Variants) > 0 {
			filtered = append(filtered, loc)
		}
	}

	// init pileup
	// pileup := make([][]*Variant, location.Len())

	// for i := range location.Len() {
	// 	pileup[i] = []*Variant{}
	// }

	// // get sorted start positions
	// starts := make([]int, 0, len(pileupMap))

	// for start := range pileupMap {
	// 	starts = append(starts, start)
	// }

	// slices.Sort(starts)

	// // assemble pileups on each start location
	// for _, start := range starts {
	// 	// sort variant types
	// 	variantTypes := make([]string, 0, len(pileupMap[start]))

	// 	for variantType := range pileupMap[start] {
	// 		variantTypes = append(variantTypes, variantType)
	// 	}

	// 	sort.Strings(variantTypes)

	// 	for _, variantType := range variantTypes {
	// 		// sort variant change
	// 		tumors := make([]string, 0, len(pileupMap[start][variantType]))

	// 		for tumor := range pileupMap[start][variantType] {
	// 			tumors = append(tumors, tumor)
	// 		}

	// 		sort.Strings(tumors)

	// 		for _, tumor := range tumors {
	// 			mutations := pileupMap[start][variantType][tumor]

	// 			for _, mutation := range mutations {
	// 				offset := start - location.Start()

	// 				pileup[offset] = append(pileup[offset], mutation)
	// 			}
	// 		}
	// 	}
	// }

	// // extract the datasets on which dataframe we are using

	// // for _, results := range search.DatasetResults {
	// // 	datasets.Add(results.Dataset)
	// // }

	return &PileupResults{Location: location, Datasets: datasets, Pileup: filtered}, nil
}

func nextYSlot(pileup *PileupLocation) (*PileupLocation, int) {
	y := pileup.Y
	pileup.Y += 1

	return pileup, y
}

func addToPileup(
	variant *Variant,
	y int,
	pileup *PileupLocation,
) *PileupLocation {

	variant.Y = y

	pileup.Variants = append(pileup.Variants, variant)

	// _, ok = (*pileupMap)[start]

	// if !ok {
	// 	(*pileupMap)[start] = make(map[string]map[string][]*Variant)
	// }

	// _, ok = (*pileupMap)[start][variant.Type]

	// if !ok {
	// 	(*pileupMap)[start][variant.Type] = make(map[string][]*Variant)
	// }

	// _, ok = (*pileupMap)[start][variant.Type][variant.Tum]

	// if !ok {
	// 	(*pileupMap)[start][variant.Type][variant.Tum] = make([]*Variant, 0, 100)
	// }

	// (*pileupMap)[start][variant.Type][variant.Tum] = append((*pileupMap)[start][variant.Type][variant.Tum], variant)

	return pileup
}

// func rowsToMutations(rows *sql.Rows) ([]*Mutation, error) {

// 	mutations := make([]*Mutation, 0, 100)

// 	defer rows.Close()

// 	for rows.Next() {
// 		var mutation Mutation

// 		err := rows.Scan(
// 			&mutation.Sample,
// 			&mutation.Chr,
// 			&mutation.Start,
// 			&mutation.End,
// 			&mutation.Ref,
// 			&mutation.Tum,
// 			&mutation.Alt,
// 			&mutation.Depth,
// 			&mutation.Type,
// 			&mutation.Vaf,
// 		)

// 		if err != nil {
// 			fmt.Println(err)
// 		}

// 		mutations = append(mutations, &mutation)
// 	}

// 	log.Debug().Msgf("all the muts %d", len(mutations))

// 	return mutations, nil
// }
