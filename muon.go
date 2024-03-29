package tmvl

import (
	"crypto/sha512"
	"math"
	"math/rand"

	"github.com/go-hep/fmom"
	"github.com/sbinet/tmvl/pumas"
)

const NMax = 3

var (
	media = []pumas.Medium{
		{
			Material:   pumas.DryAir,
			Magnetized: false,
			Density:    1.204,
			Magnet:     fmom.Vec3{0, 0, 0},
		},
		{
			Material:   pumas.StandardRock,
			Magnetized: false,
			Density:    1.5e3,
			Magnet:     fmom.Vec3{0, 0, 0},
		},
	}
)

type Muon struct {
	Energy    float64
	Direction fmom.Vec3
	Position  fmom.Vec3
	Charge    int
	Distance  float64
	Time      float64

	src rand.Source
	rng *rand.Rand
	ctx *pumas.Context
}

func NewMuon(src rand.Source) *Muon {
	mu := &Muon{
		src: src,
	}
	if mu.src == nil {
		mu.src = rand.NewSource(0)
	}
	mu.rng = rand.New(mu.src)
	return mu
}

func (mu *Muon) Run(geo Geometry) error {
	var err error

	err = mu.generate()
	if err != nil {
		return err
	}

	err = mu.propagate(geo)
	if err != nil {
		return err
	}

	return err
}

func (mu *Muon) generate() error {
	var err error
	var (
		dx   float64
		dz   float64
		norm float64
	)

	// initial energy
	mu.Energy = 5 + mu.rng.Float64()*1e4

	if mu.rng.Float64() < 0.5 {
		mu.Charge = -1
	} else {
		mu.Charge = 1
	}

	mu.Position[0] = mu.rng.Float64()*20 - 10
	mu.Position[1] = 4000
	mu.Position[2] = mu.rng.Float64()*20 - 10

	dx = mu.rng.Float64()*2 - 1
	dz = mu.rng.Float64()*2 - 1

	mu.Direction[0] = dx - mu.Position[0]
	mu.Direction[1] = 0 - mu.Position[1]
	mu.Direction[2] = dz - mu.Position[2]

	norm = math.Sqrt(mu.Direction[0]*mu.Direction[0] + mu.Direction[1]*mu.Direction[1] + mu.Direction[2]*mu.Direction[2])
	norm = 1 / norm
	mu.Direction[0] *= norm
	mu.Direction[1] *= norm
	mu.Direction[2] *= norm

	mu.ctx = pumas.New()
	mu.ctx.Rand = mu.rng
	mu.ctx.Data = mu
	mu.ctx.Step.Max = 10.0
	mu.ctx.SetLocator(
		func(ctx *pumas.Context, pos fmom.Vec3) int {
			var geo Geometry
			return geo.Medium(pos)
		},
		media,
	)

	return err
}

func (mu *Muon) propagate(geo Geometry) error {
	state, err := mu.ctx.Propagate(
		float64(mu.Charge),
		pumas.State{
			Kinetic:   mu.Energy,
			Distance:  mu.Distance,
			Time:      mu.Time,
			Position:  mu.Position,
			Direction: mu.Direction,
		},
		1, 1,
	)
	if err != nil {
		return err
	}
	mu.Energy = state.Kinetic
	mu.Distance = state.Distance
	mu.Time = state.Time
	mu.Position = state.Position
	mu.Direction = state.Direction
	return err
}

type Result struct {
	ID        [sha512.Size384]byte
	Position  fmom.Vec3
	Direction fmom.Vec3
	Energy    float64
	Time      float64
	Distance  float64
}
