//+build ignore

package main

import (
	"encoding/binary"
	"flag"
	"log"
	"math/rand"
	"os"
)

var (
	fname = flag.String("o", "prng.data", "output data file")
	nevts = flag.Int("n", 1000000, "number of events to generate")
	seed  = flag.Int64("seed", 1234, "initial seed")
)

func main() {
	flag.Parse()

	f, err := os.Create(*fname)
	if err != nil {
		log.Fatalf("error creating output file [%s]: %v\n", *fname, err)
	}
	defer f.Close()
	src := rand.NewSource(*seed)
	rng := rand.New(src)
	var buf [8]byte
	for i := 0; i < *nevts; i++ {
		v := rng.Int63()
		binary.LittleEndian.PutUint64(buf[:], uint64(v))
		_, err = f.Write(buf[:])
		if err != nil {
			log.Fatalf("error writing buffer %q (value=%v): %v\n",
				string(buf[:]),
				v,
				err,
			)
		}
	}
}
