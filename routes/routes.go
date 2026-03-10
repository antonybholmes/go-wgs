package routes

import (
	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-wgs/wgsdb"

	"github.com/antonybholmes/go-web"
	"github.com/antonybholmes/go-web/auth/token"
	"github.com/antonybholmes/go-web/middleware"
	"github.com/antonybholmes/go-wgs"

	"github.com/gin-gonic/gin"
)

type (
	MutationParams struct {
		Locations []*dna.Location
		Datasets  []string
	}

	ReqMutationParams struct {
		Locations []string `json:"locations"`
		Datasets  []string `json:"datasets"`
	}

	PileupResp struct {
		Location *dna.Location `json:"location"`
		//Info      []*mutations.MutationDBInfo `json:"info"`
		Samples   int              `json:"samples"`
		Mutations [][]*wgs.Variant `json:"mutations"`
	}

	MafResp struct {
		Location *dna.Location `json:"location"`
		//Info      []*mutations.MutationDB `json:"info"`
		Samples   int        `json:"samples"`
		Mutations [][]string `json:"mutations"`
	}
)

func ParseParamsFromPost(c *gin.Context) (*MutationParams, error) {

	var locs ReqMutationParams

	err := c.Bind(&locs)

	if err != nil {
		return nil, err
	}

	locations, err := dna.ParseLocations(locs.Locations)

	if err != nil {
		return nil, err
	}

	return &MutationParams{locations, locs.Datasets}, nil
}

func VariantDatasetsRoute(c *gin.Context) {
	middleware.JwtUserWithPermissionsRoute(c, func(c *gin.Context, isAdmin bool, user *token.AuthUserJwtClaims) {
		assembly := c.Param("assembly")

		datasets, err := wgsdb.Datasets(assembly, isAdmin, user.Permissions)

		if err != nil {
			c.Error(err)
			return
		}

		web.MakeDataResp(c, "", datasets)
	})
}

func PileupRoute(c *gin.Context) {
	searchRoute(c, "pileup")
}

func VariantsRoute(c *gin.Context) {
	searchRoute(c, "variants")
}

func searchRoute(c *gin.Context, mode string) {
	middleware.JwtUserWithPermissionsRoute(c, func(c *gin.Context, isAdmin bool, user *token.AuthUserJwtClaims) {
		assembly := c.Query("assembly")

		params, err := ParseParamsFromPost(c)

		if err != nil {
			c.Error(err)
			return
		}

		location := params.Locations[0]

		variants, err := wgsdb.Search(assembly,
			location,
			params.Datasets,
			isAdmin,
			user.Permissions)

		if err != nil {
			c.Error(err)
			return
		}

		if mode == "variants" {
			results := wgs.VariantSearchResults{
				Location: location,
				Variants: variants,
				Datasets: params.Datasets,
				//DatasetResults: make([]*DatasetResults, 0, 100)
			}

			web.MakeDataResp(c, "", results)
			return
		}

		pileup, err := wgs.GetPileup(location, params.Datasets, variants)

		if err != nil {
			c.Error(err)
			return
		}

		web.MakeDataResp(c, "", pileup)
	})
}

// func MafRoute(c *gin.Context) {
// 	 NewValidator(c).Success(func(validator *Validator) {
// 		assembly := c.Param("assembly")

// 		params, err := ParseParamsFromPost(c)

// 		if err != nil {
// 			return web.ErrorReq(err)
// 		}

// 		location := params.Locations[0]

// 		//assembly := c.Param("assembly")
// 		//name := c.Param("name")

// 		ret := MafResp{Location: location,
// 			//Info:      make([]*mutations.MutationDBInfo, len(params.Databases)),
// 			Mutations: make([][]string, location.Len())}

// 		for i := range location.Len() {
// 			ret.Mutations[i] = make([]string, 0, 10)
// 		}

// 		sampleMap := make([]map[int]struct{}, location.Len())

// 		for i := range location.Len() {
// 			sampleMap[i] = make(map[int]struct{})
// 		}

// 		for _, id := range params.Datasets {
// 			dataset, err := mutationdbcache.GetDataset(assembly, id)

// 			if err != nil {
// 				return web.ErrorReq(err)
// 			}

// 			results, err := dataset.Search(location)

// 			if err != nil {
// 				return web.ErrorReq(err)
// 			}

// 			// sum the total number of samples involved
// 			ret.Samples += len(dataset.Samples)

// 			for _, mutation := range results.Mutations {
// 				offset := mutation.Start - location.Start
// 				sample := mutation.Sample

// 				_, ok := sampleMap[offset][sample]

// 				if !ok {
// 					sampleMap[offset][sample] = struct{}{}
// 				}

// 			}

// 			//ret.Info[dbi] = db.Info
// 		}

// 		// sort each pileup
// 		for ci := range location.Len() {
// 			if len(sampleMap[ci]) > 0 {
// 				samples := make([]int, 0, len(sampleMap[ci]))

// 				for sample := range sampleMap[ci] {
// 					samples = append(samples, sample)
// 				}

// 				// sort the samples for ease of use
// 				sort.Ints(samples)

// 				ret.Mutations[ci] = append(ret.Mutations[ci], samples...)
// 			}

// 		}

// 		web.MakeDataResp(c, "", ret)
// 	})

// 	//web.MakeDataResp(c, "", mutationdbcache.GetInstance().List())
// }
