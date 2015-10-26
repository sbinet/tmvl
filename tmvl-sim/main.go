package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"math/rand"
	"os"

	"github.com/sbinet/tmvl"
)

var (
	nevts = flag.Int("nevts", 100, "number of events to compute")
	ievt  = flag.Int("ievt", 0, "first event of the event loop")
	fname = flag.String("o", "hpcsim.out", "path to output file to write")
)

type Result struct {
	Energy   float64
	Distance float64
	Time     float64
}

func main() {
	flag.Parse()

	f, err := os.Create(*fname)
	if err != nil {
		log.Fatalf("error creating output file [%s]: %v\n",
			*fname,
			err,
		)
	}

	errc := make(chan error)
	results := make(chan Result)

	go writeResults(f, results, errc)

	for i := 0; i < *nevts; i++ {
		go func(i int) {
			mu := tmvl.NewMuon(rand.NewSource(int64(i)))
			err := mu.Run(tmvl.Geometry{})
			if err != nil {
				errc <- fmt.Errorf("error propagating muon #%d: %v\n", i, err)
				return
			}
			results <- Result{
				Energy:   mu.Energy,
				Distance: mu.Distance,
				Time:     mu.Time,
			}
		}(i)
	}

	for i := 0; i < *nevts; i++ {
		err := <-errc
		if err != nil {
			log.Fatalf("error: %v\n", err)
		}
	}

}

func writeResults(f io.Writer, input <-chan Result, errc chan<- error) {
	buf := make([]byte, 3*8)
	for res := range input {
		binary.LittleEndian.PutUint64(buf[0*8:1*8], math.Float64bits(res.Energy))
		binary.LittleEndian.PutUint64(buf[1*8:2*8], math.Float64bits(res.Distance))
		binary.LittleEndian.PutUint64(buf[2*8:3*8], math.Float64bits(res.Time))
		_, err := f.Write(buf)
		if err != nil {
			err = fmt.Errorf(
				"error writing result=%#v: %v",
				res, err,
			)
		}
		errc <- err
	}
}
