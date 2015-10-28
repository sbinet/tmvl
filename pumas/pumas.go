package pumas

import (
	"errors"
	"math/rand"

	"github.com/go-hep/fmom"
)

// MaterialKind describes the kind of medium through which a particle propagates.
type MaterialKind int

const (
	StandardRock MaterialKind = iota
	DryAir
	Concrete
	Iron
	Lead
)

type PropertyKind int

const (
	Kinetic PropertyKind = iota
	Grammage
	WeightedTime
	EnergyLoss
	LarmorIntegral0
	LarmorIntegral1
	LarmorIntegral2
	LarmorIntegral3
	LarmorIntegral4
	LarmorIntegral5
	LarmorIntegral6
	LarmorIntegral7
	LarmorPhase
)

type CutKind int

const (
	MinKinetic CutKind = iota
	MaxKinetic
	MinStep
	MaxStep
	TransportRatio
	EnergyLossRatio
	MagneticRatio
)

var (
	ErrConfig            = errors.New("pumas: configuration error")
	ErrIndex             = errors.New("pumas: index error")
	ErrIO                = errors.New("pumas: I/O error")
	ErrValue             = errors.New("pumas: value error")
	ErrIncompleteMessage = errors.New("pumas: incomplete message")
)

// mode describes the propagation mode
type mode int

const (
	noMedium mode = iota
	singleMedium
	trackWithLocator
)

// Medium is a medium descriptor.
type Medium struct {
	Material   MaterialKind
	Magnetized bool
	Density    float64
	Magnet     fmom.Vec3
}

// State describes the current state of a particle propagation
type State struct {
	Medium    int
	Kinetic   float64
	Distance  float64
	Time      float64
	Position  fmom.Vec3
	Direction fmom.Vec3
}

type Recorder struct {
	First   *State
	Last    *State
	Len     int
	Scaling int
	Record  int
}

type Locator func(ctx *Context, x, y, z float64) error

type Context struct {
	mode     mode
	Medium   *Medium
	Locator  Locator
	Rand     *rand.Rand
	Recorder *Recorder
	Kinetic  struct {
		Min float64
		Max float64
	}
	Step struct {
		Min float64
		Max float64
	}
	TransportRatio  float64
	EnergyLossRatio float64
	MagneticRatio   float64

	Data interface{}
}

func (ctx *Context) Initialize() error {
	var err error
	return err
}

func (ctx *Context) Finalize() error {
	var err error
	return err
}

func (ctx *Context) Propagate(charge float64, state State, scattering, fwd int) (State, error) {
	var err error
	return state, err
}

func (ctx *Context) PropagateInMedium(medium *Medium, charge float64, state State, scattering, fwd int) (State, error) {
	var err error
	return state, err
}

func (ctx *Context) PropagateWithLocator(locator Locator, charge float64, state State, scattering, fwd int) (State, error) {
	var err error
	return state, err
}
