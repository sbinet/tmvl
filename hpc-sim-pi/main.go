package main

import (
	"encoding/binary"
	"math/rand"

	"github.com/sbinet/tmvl/sim"
)

func main() {
	sim.Main(func(app *sim.App) {
		for i := 0; i < app.NumEvents(); i++ {
			go simulate(i, app.Results())
		}
	})
}

func simulate(ievt int, results chan<- sim.Result) {
	src := rand.NewSource(int64(ievt))
	rng := rand.New(src)

	total := 0
	inside := 0

	for total = 0; total < 100000; total++ {
		x := rng.Float64()
		y := rng.Float64()
		if x*x+y*y < 1.0 {
			inside++
		}
	}

	buf := make([]byte, 2*8)
	binary.LittleEndian.PutUint64(buf[:8], uint64(total))
	binary.LittleEndian.PutUint64(buf[8:], uint64(inside))
	results <- sim.Result{ID: int64(ievt), Data: buf}
}
