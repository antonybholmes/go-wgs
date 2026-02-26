package wgsdb

import (
	"sync"

	"github.com/antonybholmes/go-dna"
	wgs "github.com/antonybholmes/go-wgs"
)

var instance *wgs.WGSDB
var once sync.Once

func InitDB(dir string) (*wgs.WGSDB, error) {
	once.Do(func() {
		instance = wgs.NewWGSDB(dir)
	})

	return instance, nil
}

func GetInstance() *wgs.WGSDB {
	return instance
}

func Datasets(assembly string, isAdmin bool, permissions []string) ([]*wgs.Dataset, error) {
	return instance.Datasets(assembly, isAdmin, permissions)
}

func Dir() string {
	return instance.Dir()
}

// func Dataset(datasetId string, isAdmin bool, permissions []string) (*mutations.Dataset, error) {
// 	return instance.Dataset(datasetId, isAdmin, permissions)
// }

func Search(assembly string, location *dna.Location, publicIds []string, isAdmin bool, permissions []string) (*wgs.SearchResults, error) {
	return instance.Search(assembly, location, publicIds, isAdmin, permissions)
}
