package wgs

import (
	"database/sql"
	"fmt"
	"path/filepath"

	"strings"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-sys"
	"github.com/antonybholmes/go-web"
	"github.com/antonybholmes/go-web/auth/sqlite"
	"github.com/rs/zerolog/log"
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
		Mutations   int       `json:"mutations"`
	}

	Variant struct {
		Dataset     string  `json:"dataset"`
		Sample      string  `json:"sample"`
		Chr         string  `json:"chr"`
		Ref         string  `json:"ref"`
		Tum         string  `json:"tum"`
		GeneSymbol  string  `json:"geneSymbol"`
		Type        string  `json:"type"` // SNV, DEL, INS
		Start       int     `json:"start"`
		End         int     `json:"end"`
		TRefCount   int     `json:"tRefCount"`
		TAltCount   int     `json:"tAltCount"`
		TDepth      int     `json:"tDepth"`
		HGVSc       string  `json:"hgvsC,omitempty"`
		HGVSp       string  `json:"hgvsP,omitempty"`
		Consequence string  `json:"consequence,omitempty"`
		Vaf         float64 `json:"vaf"`
		Y           int     `json:"y"` // fast access for plotting
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
		Pos      int        `json:"pos"`
		// keep track of max y at position
		y int `json:"-"`
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

	WGSDB struct {
		db  *sql.DB
		dir string
	}
)

const (
	DatasetsSql = `SELECT DISTINCT
		d.public_id AS dataset_id,
		g.name AS genome,
		a.name AS assembly,
		ins.name AS institution,
		d.name,
		d.short_name,
		d.mutations,
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
		JOIN samples s ON s.dataset_id = d.id
		JOIN dataset_permissions dp ON d.id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE 
			<<PERMISSIONS>>
			AND LOWER(a.name) = :assembly
		ORDER BY
			d.name, 
			s.name`

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

	FindVariantsSql = `SELECT DISTINCT
		d.public_id AS dataset_id,
		s.public_id AS sample_id,
		c.name AS chr, 
		v.start, 
		v.end,
		v.ref, 
		v.tum, 
		vt.name AS variant_type,
		'' AS gene_symbol,
		sv.t_depth - sv.t_alt_count AS ref_count, 
		sv.t_alt_count, 
		sv.t_depth,
		v.hgvs_c,
		v.hgvs_p,
		v.consequence,
		sv.vaf
		FROM sample_variants sv
		JOIN variants v ON sv.variant_id = v.id
		JOIN chromosomes c ON c.id = v.chr_id
		JOIN variant_types vt ON vt.id = v.variant_type_id
		JOIN samples s ON sv.sample_id = s.id
		JOIN datasets d ON s.dataset_id = d.id
		JOIN dataset_permissions dp ON s.dataset_id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE
			<<PERMISSIONS>>
			AND <<DATASETS>>
			AND c.name = :chr AND v.end >= :start AND v.start <= :end
		ORDER BY vt.id, v.start, d.name`
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

func (mutation *Variant) Clone() *Variant {
	var ret Variant = Variant{Chr: mutation.Chr,
		Start:     mutation.Start,
		End:       mutation.End,
		Ref:       mutation.Ref,
		Tum:       mutation.Tum,
		TRefCount: mutation.TRefCount,
		TAltCount: mutation.TAltCount,
		TDepth:    mutation.TDepth,
		Type:      mutation.Type,
		Vaf:       mutation.Vaf,
		Sample:    mutation.Sample,
		Dataset:   mutation.Dataset,
		Y:         mutation.Y,
	}

	return &ret
}

func NewWGSDB(dir string) *WGSDB {

	db := sys.Must(sql.Open(sys.Sqlite3DB, filepath.Join(dir, "wgs.db"+sys.SqliteDSN)))

	return &WGSDB{dir: dir, db: db}
}

func (mdb *WGSDB) Dir() string {
	return mdb.dir
}

func (mdb *WGSDB) Datasets(assembly string, isAdmin bool, permissions []string) ([]*Dataset, error) {
	namedArgs := []any{sql.Named("assembly", web.FormatParam(assembly))}

	query := sqlite.MakePermissionsSql(DatasetsSql, isAdmin, permissions, &namedArgs)

	rows, err := mdb.db.Query(query, namedArgs...)

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

		err := rows.Scan(&dataset.PublicId,
			&dataset.Genome,
			&dataset.Assembly,
			&dataset.Institution,
			&dataset.Name,
			&dataset.ShortName,
			&dataset.Mutations,
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

func (mdb *WGSDB) Search(assembly string,
	location *dna.Location,
	datasetIds []string,
	isAdmin bool,
	permissions []string) ([]*Variant, error) {

	namedArgs := []any{
		sql.Named("chr", location.Chr()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End())}

	query := sqlite.MakePermissionsSql(FindVariantsSql, isAdmin, permissions, &namedArgs)

	query = MakeInDatasetsSql(query, datasetIds, &namedArgs)

	rows, err := mdb.db.Query(query, namedArgs...)

	log.Debug().Msgf("querying mutations db with %v %s", isAdmin, query)

	if err != nil {
		log.Debug().Msgf("error querying mutations db: %s", err)
		return nil, err
	}

	defer rows.Close()

	results := make([]*Variant, 0, 100)

	//var currentDatasetResults *DatasetResults

	for rows.Next() {
		var variant Variant

		err := rows.Scan(
			&variant.Dataset,
			&variant.Sample,
			&variant.Chr,
			&variant.Start,
			&variant.End,
			&variant.Ref,
			&variant.Tum,
			&variant.Type,
			&variant.GeneSymbol,
			&variant.TRefCount,
			&variant.TAltCount,
			&variant.TDepth,
			&variant.HGVSc,
			&variant.HGVSp,
			&variant.Consequence,
			&variant.Vaf,
		)

		//log.Debug().Msgf("found mutation %v", mutation)

		if err != nil {
			return nil, err
		}

		// since ordered by dataset we can group here
		// if currentDatasetResults == nil || currentDatasetResults.Dataset != variant.Dataset {
		// 	// create a new stores for this dataset and set to current
		// 	currentDatasetResults = &DatasetResults{Dataset: variant.Dataset, Variants: make([]*Variant, 0, 100)}
		// 	results.DatasetResults = append(results.DatasetResults, currentDatasetResults)
		// }

		results = append(results, &variant)
	}

	return results, nil
}

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
		pileups[i] = &PileupLocation{Start: location.Start() + i, Pos: i, Variants: make([]*Variant, 0, 100)}
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
					pileups[pos].y = pileups[pos-1].y
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
					pileups[pos].y = pileups[pos-1].y
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
	y := pileup.y
	pileup.y += 1

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
