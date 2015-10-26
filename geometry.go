package tmvl

type Geometry struct {
}

func (geo *Geometry) Medium(x, y, z float64) int {
	// FIXME(sbinet): devise a more realistic implementation
	// so far: a big 3d-block of 5m x 2000m x 5m
	//         centered around the y-axis
	//         density of 1.66 g/cm^3
	// otherwise: open sky density of 0.001204 g/cm^3
	// if y < 0, we reached the plane of the detector.
	switch {
	case y > 1.0 && y < 2001 &&
		x > -2.5 && x < 2.5 &&
		z > -2.5 && z < 2.5:
		return 1
	case y < 0.0:
		return -1
	}
	return 0
}
