package pumas

import (
	"math/rand"
)

type Context struct {
	Rand *rand.Rand
	Data interface{}
}

func (ctx *Context) Propagate(charge, ene float64, pos, dir [3]float64, dist, t float64, v0, v1 int) (float64, float64, float64) {

	return ene, dist, t
}
