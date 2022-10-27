//+build ignore

package main

import "flag"

var (
	nevts  = flag.Int("nevts", 100, "number of events to compute")
	ievt   = flag.Int("ievt", 0, "first event of the event loop")
	fname  = flag.String("o", "hpcsim.out", "path to output file to write")
	simpkg = flag.String("sim", "", "fully qualified name of the package to run")
)

type SimContext struct {
	NEvts int
	IEvt  int
	Data  interface{} // user data
}

type Sim interface {
	SimStart(ctx *SimContext) error
	SimStop(ctx *SimContext) error
}

type Runner interface {
	RunStart(ctx *SimContext) error
	RunStop(ctx *SimContext) error
}

type Context struct{}

type Eventer interface {
	EventStart(sim *SimContext, ctx *Context) error
	EventRun(sim *SimContext, ctx *Context) error
	EventStop(sim *SimContext, ctx *Context) error
}

type Result struct {
	ID   [48]byte
	Data []byte
}

func main() {
	flag.Parse()

}
