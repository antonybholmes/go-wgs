package wgsdb

import (
	"sync"

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

func Instance() *wgs.WGSDB {
	return instance
}
