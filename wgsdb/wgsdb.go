package wgsdb

import (
	"sync"

	"github.com/antonybholmes/go-dna"
	"github.com/antonybholmes/go-web/auth/token"
	wgs "github.com/antonybholmes/go-wgs"
)

var instance *wgs.WGSDB
var once sync.Once

func InitDB(file string) (*wgs.WGSDB, error) {
	once.Do(func() {
		instance = wgs.NewWGSDB(file)
	})

	return instance, nil
}

func GetInstance() *wgs.WGSDB {
	return instance
}

func Datasets(assembly string, isAdmin bool, permissions []string) ([]*wgs.Dataset, error) {
	return instance.Datasets(assembly, isAdmin, permissions)
}

// func Dir() string {
// 	return instance.Dir()
// }

// func Dataset(datasetId string, isAdmin bool, permissions []string) (*mutations.Dataset, error) {
// 	return instance.Dataset(datasetId, isAdmin, permissions)
// }

func Search(assembly string, location *dna.Location, datasetIds []string, isAdmin bool, permissions []string) ([]*wgs.Variant, error) {
	return instance.Search(assembly, location, datasetIds, isAdmin, permissions)
}

func MAF(assembly string, location *dna.Location, datasetIds []string, isAdmin bool, user *token.AuthUserJwtClaims) (*wgs.MAFResults, error) {
	return instance.MAF(assembly, location, datasetIds, isAdmin, user)
}
