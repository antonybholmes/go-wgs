package mutations

import (
	"database/sql"
	"fmt"
	"path/filepath"

	"slices"
	"sort"
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

	Mutation struct {
		Dataset string  `json:"dataset"`
		Sample  string  `json:"sample"`
		Chr     string  `json:"chr"`
		Ref     string  `json:"ref"`
		Tum     string  `json:"tum"`
		Type    string  `json:"type"`
		Start   int     `json:"start"`
		End     int     `json:"end"`
		Alt     int     `json:"tAltCount"`
		Depth   int     `json:"tDepth"`
		Vaf     float64 `json:"vaf"`
	}

	DatasetResults struct {
		Dataset   string      `json:"dataset"`
		Mutations []*Mutation `json:"mutations"`
	}

	PileupResults struct {
		Location *dna.Location `json:"location"`
		Datasets []string      `json:"datasets"`
		Pileup   [][]*Mutation `json:"pileup"`
	}

	SearchResults struct {
		Location       *dna.Location     `json:"location"`
		DatasetResults []*DatasetResults `json:"results"`
	}

	MutationsDB struct {
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

	FindMutationsSql = `SELECT
		d.public_id AS dataset_id,
		s.public_id AS sample_id,
		c.name, 
		m.start, 
		m.end, 
		m.ref, 
		m.tum, 
		m.t_alt_count, 
		m.t_depth, 
		m.variant_type,
		m.vaf
		FROM mutations m
		JOIN chromosomes c ON c.id = m.chr_id
		JOIN samples s ON m.sample_id = s.id
		JOIN datasets d ON s.dataset_id = d.id
		JOIN dataset_permissions dp ON s.dataset_id = dp.dataset_id
		JOIN permissions p ON dp.permission_id = p.id
		WHERE
			<<PERMISSIONS>>
			AND <<DATASETS>>
			AND c.name = :chr AND m.start >= :start AND  m.end <= :end
		ORDER BY d.public_id, c.name, m.start, m.end, m.variant_type`
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

func (mutation *Mutation) Clone() *Mutation {
	var ret Mutation = Mutation{Chr: mutation.Chr,
		Start:  mutation.Start,
		End:    mutation.End,
		Ref:    mutation.Ref,
		Tum:    mutation.Tum,
		Alt:    mutation.Alt,
		Depth:  mutation.Depth,
		Type:   mutation.Type,
		Vaf:    mutation.Vaf,
		Sample: mutation.Sample,
	}

	return &ret
}

func NewMutationsDB(dir string) *MutationsDB {

	db := sys.Must(sql.Open(sys.Sqlite3DB, filepath.Join(dir, "mutations.db"+sys.SqliteReadOnlySuffix)))

	return &MutationsDB{dir: dir, db: db}
}

func (mdb *MutationsDB) Dir() string {
	return mdb.dir
}

func (mdb *MutationsDB) Datasets(assembly string, isAdmin bool, permissions []string) ([]*Dataset, error) {
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

func (mdb *MutationsDB) Search(assembly string,
	location *dna.Location,
	datasetIds []string,
	isAdmin bool,
	permissions []string) (*SearchResults, error) {
	results := SearchResults{Location: location, DatasetResults: make([]*DatasetResults, 0, len(datasetIds))}

	namedArgs := []any{
		sql.Named("chr", location.Chr()),
		sql.Named("start", location.Start()),
		sql.Named("end", location.End())}

	query := sqlite.MakePermissionsSql(FindMutationsSql, isAdmin, permissions, &namedArgs)

	query = MakeInDatasetsSql(query, datasetIds, &namedArgs)

	rows, err := mdb.db.Query(query, namedArgs...)

	log.Debug().Msgf("querying mutations db with %v %s", isAdmin, query)

	if err != nil {
		log.Debug().Msgf("error querying mutations db: %s", err)
		return nil, err
	}

	defer rows.Close()

	var currentDatasetResults *DatasetResults

	for rows.Next() {
		var mutation Mutation

		err := rows.Scan(
			&mutation.Dataset,
			&mutation.Sample,
			&mutation.Chr,
			&mutation.Start,
			&mutation.End,
			&mutation.Ref,
			&mutation.Tum,
			&mutation.Alt,
			&mutation.Depth,
			&mutation.Type,
			&mutation.Vaf,
		)

		//log.Debug().Msgf("found mutation %v", mutation)

		if err != nil {
			return nil, err
		}

		// since ordered by dataset we can group here
		if currentDatasetResults == nil || currentDatasetResults.Dataset != mutation.Dataset {
			// create a new stores for this dataset and set to current
			currentDatasetResults = &DatasetResults{Dataset: mutation.Dataset, Mutations: make([]*Mutation, 0, 100)}
			results.DatasetResults = append(results.DatasetResults, currentDatasetResults)
		}

		currentDatasetResults.Mutations = append(currentDatasetResults.Mutations, &mutation)
	}

	return &results, nil
}

func GetPileup(search *SearchResults) (*PileupResults, error) {
	// first lets fix deletions and insertions
	for _, datasetResults := range search.DatasetResults {
		for _, mutation := range datasetResults.Mutations {
			// change for sorting purposes so that ins always comes last
			switch mutation.Type {
			case "INS":
				mutation.Type = "2:INS"
				// modify the output so that is begins with a caret to indicate
				// an insertion
				mutation.Tum = fmt.Sprintf("^%s", mutation.Tum)
			case "DEL":
				mutation.Type = "3:DEL"
			default:
				mutation.Type = "1:SNP"
			}
		}
	}

	// put together by position, type, tum

	pileupMap := make(map[int]map[string]map[string][]*Mutation)

	for _, datasetResults := range search.DatasetResults {
		for _, mutation := range datasetResults.Mutations {

			mutation.Dataset = datasetResults.Dataset

			switch mutation.Type {
			case "3:DEL":
				for i := range mutation.End - mutation.Start + 1 {
					addToPileupMap(mutation.Start+i, mutation, &pileupMap)
				}
			case "2:INS":
				addToPileupMap(mutation.Start, mutation, &pileupMap)
			default:
				// deal with concatenated snps e.g. ACG>TGT
				//tum := []rune(mutation.Tum)
				for i, c := range mutation.Tum {
					// clone and change tumor
					mut2 := mutation.Clone()
					mut2.Tum = string(c)
					addToPileupMap(mut2.Start+i, mut2, &pileupMap)
				}
			}
		}
	}

	location := search.Location

	// init pileup
	pileup := make([][]*Mutation, location.Len())

	for i := range location.Len() {
		pileup[i] = []*Mutation{}
	}

	// get sorted start positions
	starts := make([]int, 0, len(pileupMap))

	for start := range pileupMap {
		starts = append(starts, start)
	}

	slices.Sort(starts)

	// assemble pileups on each start location
	for _, start := range starts {
		// sort variant types
		variantTypes := make([]string, 0, len(pileupMap[start]))

		for variantType := range pileupMap[start] {
			variantTypes = append(variantTypes, variantType)
		}

		sort.Strings(variantTypes)

		for _, variantType := range variantTypes {
			// sort variant change
			tumors := make([]string, 0, len(pileupMap[start][variantType]))

			for tumor := range pileupMap[start][variantType] {
				tumors = append(tumors, tumor)
			}

			sort.Strings(tumors)

			for _, tumor := range tumors {
				mutations := pileupMap[start][variantType][tumor]

				for _, mutation := range mutations {
					offset := start - location.Start()

					pileup[offset] = append(pileup[offset], mutation)
				}
			}
		}
	}

	// extract the datasets on which dataframe we are using
	datasets := make([]string, 0, len(search.DatasetResults))

	for _, results := range search.DatasetResults {
		datasets = append(datasets, results.Dataset)
	}

	return &PileupResults{Location: location, Datasets: datasets, Pileup: pileup}, nil
}

func addToPileupMap(start int,
	mutation *Mutation,
	pileupMap *map[int]map[string]map[string][]*Mutation) {

	_, ok := (*pileupMap)[start]

	if !ok {
		(*pileupMap)[start] = make(map[string]map[string][]*Mutation)
	}

	_, ok = (*pileupMap)[start][mutation.Type]

	if !ok {
		(*pileupMap)[start][mutation.Type] = make(map[string][]*Mutation)
	}

	_, ok = (*pileupMap)[start][mutation.Type][mutation.Tum]

	if !ok {
		(*pileupMap)[start][mutation.Type][mutation.Tum] = make([]*Mutation, 0, 100)
	}

	(*pileupMap)[start][mutation.Type][mutation.Tum] = append((*pileupMap)[start][mutation.Type][mutation.Tum], mutation)
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
