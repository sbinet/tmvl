package main

import (
	"bufio"
	"encoding/binary"
	"flag"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"runtime/pprof"
	"runtime/trace"
)

var (
	nthreads = flag.Int("t", 1, "number of concurrent goroutines")
	nevts    = flag.Int("nevts", 100, "number of events to compute")
	ievt     = flag.Int("ievt", 0, "first event of the event loop")
	fname    = flag.String("o", "hpcsim.out", "path to output file to write")
	simpkg   = flag.String("sim", "", "fully qualified name of the package to run")

	doprof  = flag.Bool("prof", false, "enable CPU profiling")
	dotrace = flag.Bool("trace", false, "enable tracing")
)

type Result struct {
	ID   [48]byte
	Data []byte
}

type Context struct {
	Total  float64
	Inside float64
}

type Sim struct {
	NEvts int
	IEvt  int
}

type App struct {
	nevts int
	ievt  int
	fname string
	procs []Proc

	errc chan error
	resc chan Result
}

type Proc interface {
}

func main() {
	flag.Parse()

	if *doprof {
		fprof, err := os.Create("prof.out")
		if err != nil {
			log.Fatalf("error creating pprof output file: %v\n", err)
		}
		defer fprof.Close()
		pprof.StartCPUProfile(fprof)
		defer pprof.StopCPUProfile()
	}

	if *dotrace {
		ftrace, err := os.Create("trace.out")
		if err != nil {
			log.Fatalf("error creating trace output file: %v\n", err)
		}
		defer ftrace.Close()
		err = trace.Start(ftrace)
		if err != nil {
			log.Fatalf("error starting tracer: %v\n", err)
		}
		defer trace.Stop()
	}

	sim := Sim{
		NEvts: *nevts,
		IEvt:  *ievt,
	}

	err := sim.Initialize()
	if err != nil {
		log.Fatalf("error initializing simulation: %v\n", err)
	}

	f, err := os.Create(*fname)
	if err != nil {
		log.Fatalf("error creating output file [%s]: %v\n",
			*fname,
			err,
		)
	}
	defer f.Close()

	errc := make(chan error, *nthreads)
	results := make(chan Result, *nthreads)
	go writeResults(f, results, errc)

	//go sim.Run(results)
	for i := 0; i < sim.NEvts; i++ {
		go sim.Simulate(i, results)
	}

	for i := 0; i < sim.NEvts; i++ {
		err := <-errc
		if err != nil {
			log.Fatalf("error: %v\n", err)
		}
	}
}

func (sim *Sim) Initialize() error {
	var err error
	return err
}

func (sim *Sim) Simulate(ievt int, results chan<- Result) {
	src := rand.NewSource(int64(ievt))
	rng := rand.New(src)

	total := 0
	inside := 0

	for total = 0; total < 10000; total++ {
		x := rng.Float64()
		y := rng.Float64()
		if x*x+y*y < 1.0 {
			inside++
		}
	}

	buf := make([]byte, 2*8)
	binary.LittleEndian.PutUint64(buf[:8], uint64(total))
	binary.LittleEndian.PutUint64(buf[8:], uint64(inside))
	results <- Result{Data: buf}
}

func writeResults(f io.Writer, input <-chan Result, errc chan<- error) {
	w := bufio.NewWriter(f)
	defer w.Flush()
	for res := range input {
		_, err := f.Write(res.Data)
		if err != nil {
			err = fmt.Errorf(
				"error writing result-id=%q: %v",
				string(res.ID[:]), err,
			)
		}
		errc <- err
	}
}
