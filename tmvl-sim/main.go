package main

import (
	"encoding/binary"
	"fmt"
	"math"
	"math/rand"

	"github.com/sbinet/tmvl"
	"github.com/sbinet/tmvl/sim"
)

func main() {
	sim.Main(func(app *sim.App) {
		for i := 0; i < app.NumEvents(); i++ {
			go simulate(i, app)
		}
	})
}

func simulate(i int, app *sim.App) {
	mu := tmvl.NewMuon(rand.NewSource(int64(i)))
	err := mu.Run(tmvl.Geometry{})
	if err != nil {
		app.Errors() <- fmt.Errorf("error propagating muon #%d: %v\n", i, err)
		return
	}
	buf := make([]byte, 3*8)
	binary.LittleEndian.PutUint64(buf[0*8:1*8], math.Float64bits(mu.Energy))
	binary.LittleEndian.PutUint64(buf[1*8:2*8], math.Float64bits(mu.Distance))
	binary.LittleEndian.PutUint64(buf[2*8:3*8], math.Float64bits(mu.Time))
	app.Results() <- sim.Result{
		ID:   int64(i),
		Data: buf,
	}
}
