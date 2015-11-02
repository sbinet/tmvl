package pumas

import (
	"errors"
	"fmt"
	"log"
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

type Locator func(ctx *Context, pos fmom.Vec3) int

type Context struct {
	mode     mode
	media    []Medium
	locator  Locator
	Rand     *rand.Rand
	Recorder *Recorder
	Kinetic  struct {
		Min float64 // in GeV
		Max float64 // in GeV
	}
	Step struct {
		Min float64 // in meters
		Max float64 // in meters
	}
	TransportRatio  float64
	EnergyLossRatio float64
	MagneticRatio   float64

	Data interface{}
}

func New() *Context {
	ctx := &Context{
		mode:     noMedium,
		media:    nil,
		locator:  nil,
		Rand:     nil, // FIXME(sbinet): use the default from math/rand
		Recorder: nil,
	}
	ctx.Step.Min = 1e-3 // meters
	ctx.Step.Max = 10   // meters
	ctx.TransportRatio = 1e-3
	ctx.EnergyLossRatio = 1e-3
	ctx.MagneticRatio = 1e-3
	return ctx
}

func (ctx *Context) Initialize() error {
	var err error
	return err
}

func (ctx *Context) Finalize() error {
	var err error
	return err
}

// Propagate propagates through a pre-configured set of media using either:
//  - propagation through a uniform medium
//  - propagation through a set of media described by a locator function.
func (ctx *Context) Propagate(charge float64, state State, scattering, fwd int) (State, error) {
	var err error
	switch ctx.mode {
	case noMedium:
		return state, fmt.Errorf("pumas: no propagation medium")
	case singleMedium:
		return ctx.propagateInMedium(ctx.media, charge, state, scattering, fwd)
	case trackWithLocator:
		return ctx.propagateWithLocator(ctx.locator, ctx.media, charge, state, scattering, fwd)
	default:
		return state, fmt.Errorf(
			"pumas: invalid propagation mode (mode=%d)",
			ctx.mode,
		)
	}
	return state, err
}

func (ctx *Context) propagateInMedium(media []Medium, charge float64, state State, scattering, fwd int) (State, error) {
	var err error
	panic("not implemented")
	return state, err
}

// propagateWithLocator propagates through a set of media described by a locator
// function.
//
// Reverse propagation allows to compute the required minimum kinetic energy for
// reaching the starting point from the world boundary.
// The corresponding distance and proper time are also given.
// Note: time is negative for a reverse propagation.
func (ctx *Context) propagateWithLocator(locator Locator, media []Medium, charge float64, state State, scattering, fwd int) (State, error) {
	var err error

	// check initial position
	idx := locator(ctx, state.Position)
	if idx < 0 {
		return state, err
	}

	// FIXME(sbinet): implement track recorder
	// if record { ... }

	// step through the media
	sign := 2 * float64(fwd-1)
	medium := media[idx]
	dist := 0.0
	xi := grammage(medium.Material, state.Kinetic)

	// check kinetic cuts
	length, ok := ctx.lengthLimit(fwd, state.Kinetic, xi, medium.Material, medium.Density)
	if !ok {
		return state, err
	}

	step := 1
proploop:
	for {
		log.Printf("--- step %d...\n", step)
		x := xi - sign*dist*medium.Density
		k := 0.0
		if x > 0 {
			k = kinetic(medium.Material, x)
		} else {
			k = -1
		}

		if k <= 0 {
			state.Kinetic = 0
			state.Time += weightedTime(medium.Material, state.Kinetic) / medium.Density * sign
			state.Distance += grammage(medium.Material, state.Kinetic) / medium.Density
			break proploop
		}

		istate := state
		istate.Kinetic = k
		istate.Distance = dist
		newIdx := ctx.stepThrough(locator, media, length, idx, charge, &istate, scattering, fwd)
		state.Position = istate.Position
		state.Direction = istate.Direction
		if newIdx != idx {
			state.Distance += istate.Distance
			X := xi - sign*istate.Distance*medium.Density
			if X > 0 {
				state.Kinetic = kinetic(medium.Material, X)
			} else {
				state.Kinetic = 0
			}
			state.Time += (weightedTime(medium.Material, k) - weightedTime(medium.Material, state.Kinetic)) / medium.Density
			idx = newIdx
			if idx < 0 {
				break proploop
			}

			// FIXME(sbinet): implement track recorder
			// if record { ... }

			medium = media[idx]
			dist = 0.0
			xi = grammage(medium.Material, state.Kinetic)
			length, ok = ctx.lengthLimit(fwd, state.Kinetic, xi, medium.Material, medium.Density)
			if !ok {
				break proploop
			}
		} else if ctx.Recorder != nil { // FIXME(sbinet): implement recording
		}
		step++
	}

	// FIXME(sbinet): implement recording
	// if recorder { ... }

	return state, err
}

func (ctx *Context) SetLocator(locator Locator, media []Medium) {
	ctx.media = media
	ctx.locator = locator
	ctx.mode = trackWithLocator
}

func (ctx *Context) SetMedium(media []Medium) {
	ctx.media = media
	ctx.locator = nil
	ctx.mode = singleMedium
}
