package sim

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime/pprof"
	"runtime/trace"
)

var (
	app = &App{}
)

func Main(f func(*App)) {
	main(f)
}

type Result struct {
	ID   int64
	Data []byte
}

type App struct {
	nconc int // number of concurrent goroutines
	nevts int
	ievt  int
	fname string

	fprof  string
	ftrace string

	errc chan error
	resc chan Result
}

func main(f func(*App)) {
	flag.Parse()

	if app.fprof != "" {
		fprof, err := os.Create(app.fprof)
		if err != nil {
			log.Fatalf("error creating pprof output file [%s]: %v\n",
				app.fprof,
				err,
			)
		}
		defer fprof.Close()
		pprof.StartCPUProfile(fprof)
		defer pprof.StopCPUProfile()
	}

	if app.ftrace != "" {
		ftrace, err := os.Create(app.ftrace)
		if err != nil {
			log.Fatalf("error creating trace output file: %v\n", err)
		}
		defer ftrace.Close()
		trace.Start(ftrace)
		defer trace.Stop()
	}

	fout, err := os.Create(app.fname)
	if err != nil {
		log.Fatalf("error creating output file [%s]: %v\n",
			app.fname,
			err,
		)
	}
	defer fout.Close()

	go app.write(fout)

	app.errc = make(chan error, app.nconc)
	app.resc = make(chan Result, app.nconc)

	go f(app)

	for i := 0; i < app.nevts; i++ {
		err := <-app.errc
		if err != nil {
			log.Fatalf("error: %v\n", err)
		}
	}
	close(app.resc)
	close(app.errc)

	err = fout.Close()
	if err != nil {
		log.Fatalf("error closing output file [%s]: %v\n",
			app.fname,
			err,
		)
	}
}

func (app *App) Errors() chan<- error {
	return app.errc
}

func (app *App) Results() chan<- Result {
	return app.resc
}

func (app *App) NumEvents() int {
	return app.nevts
}

func init() {
	flag.IntVar(&app.nconc, "nprocs", 1, "number of concurrent goroutines")
	flag.IntVar(&app.nevts, "nevts", 100, "number of events to process")
	flag.IntVar(&app.ievt, "ievt", 0, "first event of the event loop")
	flag.StringVar(&app.fname, "o", "hpcsim.out", "path to output file to store results")
	flag.StringVar(&app.fprof, "cpu-profile", "", "enable CPU profiling")
	flag.StringVar(&app.ftrace, "ftrace", "", "enable tracing")
}

func (app *App) write(f io.Writer) {
	w := bufio.NewWriter(f)
	defer w.Flush()
	for res := range app.resc {
		_, err := w.Write(res.Data)
		if err != nil {
			err = fmt.Errorf(
				"error writing result-id=%d: %v",
				res.ID, err,
			)
		}
		app.errc <- err
	}
}
