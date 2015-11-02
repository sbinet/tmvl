package pumas

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/andrebq/gas"
	"github.com/go-hep/fmom"
)

const (
	numPDG       = 146
	numMaterials = 5
	numElements  = 14
	numTaylor    = 8
	maxComps     = 10
)

type elementType struct {
	z    float64 // atomic number
	a    float64 // atomic mass in g/mol
	name string
}

type componentType struct {
	elem int     // element index
	frac float64 // mass fraction
}

// Physical constants
const (
	physLifetime float64 = 658.654     // in m
	physKLarmor  float64 = 0.299792458 // in m^-1 GeV/c T^-1
	physEscat2   float64 = 184.96e-6   // in GeV^2
)

var (
	elements [numElements]elementType
	physK    [numPDG]float64
)

type materialDef struct {
	numElem int
	comp    []componentType
	x       [numPDG]float64
	t       [numPDG]float64
	dE      [numPDG]float64
	lp      [numPDG]float64 // Larmor's phase
	li      [numTaylor][numPDG]float64
	amax    float64
	bmax    float64
	density float64
	radlen  float64
	name    string
	kmi     float64
	kcrit   float64
}

var (
	materials = [numMaterials]materialDef{}
)

var (
	muon muonDef
)

type muonDef struct {
	Mass float64
}

// loadTables loads the tabulated CSDA data into memory
func loadTables(dir string) error {
	var err error
	fnames := [numMaterials]string{
		"standard_rock.dat",
		"dry_air.dat",
		"concrete.dat",
		"iron.dat",
		"lead.dat",
	}

	materials[StandardRock].numElem = 1
	materials[DryAir].numElem = 4
	materials[Concrete].numElem = 10
	materials[Iron].numElem = 1
	materials[Lead].numElem = 1

	materials[StandardRock].radlen = 265.4
	materials[DryAir].radlen = 366.2
	materials[Concrete].radlen = 265.7
	materials[Iron].radlen = 138.4
	materials[Lead].radlen = 63.7

	elements = [numElements]elementType{
		{
			name: "Hydrogen:Gas",
			z:    1.0,
			a:    1.00795, /* g/mol            */
		},
		{
			name: "Carbon:Amorphous",
			z:    6.,
			a:    12.0108, /* g/mol            */
		},
		{
			name: "Nitrogen:Gas",
			z:    7.,
			a:    14.0067, /* g/mol            */
		},
		{
			name: "Oxygen:Gas",
			z:    8.,
			a:    15.9994, /* g/mol            */
		},
		{
			name: "Sodium",
			z:    11.,
			a:    22.9898, /* g/mol            */
		},
		{
			name: "Rock:Standard",
			z:    11.,
			a:    22., /* g/mol            */
		},
		{
			name: "Magnesium",
			z:    12.,
			a:    24.3051, /* g/mol            */
		},
		{
			name: "Aluminium",
			z:    13.,
			a:    26.9815, /* g/mol            */
		},
		{
			name: "Silicon",
			z:    14.,
			a:    28.0855, /* g/mol            */
		},
		{
			name: "Argon:Gas",
			z:    18.,
			a:    39.948, /* g/mol            */
		},
		{
			name: "Potassium",
			z:    19.,
			a:    39.0983, /* g/mol            */
		},
		{
			name: "Calcium",
			z:    20.,
			a:    40.0784, /* g/mol            */
		},
		{
			name: "Iron",
			z:    26.,
			a:    55.8452, /* g/mol            */
		},
		{
			name: "Lead",
			z:    82.,
			a:    207.2, /* g/mol            */
		},
	}

	if dir == "" {
		dir = os.Getenv("PUMAS_CSDA_DIR")
	}

	if dir == "" {
		return fmt.Errorf("pumas: no $PUMAS_CSDA_DIR defined")
	}

	for imed, fname := range fnames {
		f, err := os.Open(filepath.Join(dir, fname))
		if err != nil {
			return fmt.Errorf(
				"pumas: could not open file #%d [%s]: %v\n",
				imed,
				filepath.Join(dir, fname),
				err,
			)
		}
		defer f.Close()

		var (
			//	itmp int
			line string
		)

		buf := make([]byte, 512)
		_, err = f.Read(buf[:41])
		if err != nil {
			return err
		}
		mass := 0.0
		// muon mass value used for the tables
		_, err = fmt.Fscanf(f, "%f MeV\n", &mass)
		if err != nil {
			return err
		}
		if imed == 0 {
			muon.Mass = mass * 1e-3 // convert from MeV to GeV
		} else {
			if mass*1e-3 != muon.Mass {
				return fmt.Errorf(
					"pumas: muon mass inconsistency (got=%v, want=%v)",
					mass*1e-3,
					muon.Mass,
				)
			}
		}

		mat := &materials[imed]
		r := bufio.NewReader(f)
		line, err = r.ReadString('\n')
		if err != nil {
			return err
		}

		// medium name
		tokens := strings.Split(line, ":")
		mat.name = strings.TrimSpace(tokens[1])

		// density
		line, err = r.ReadString('\n')
		if err != nil {
			return err
		}
		mat.density, err = parseDensity(line)
		if err != nil {
			log.Fatalf("error parsing density [%s]: %v\n",
				fname,
				err,
			)
			return err
		}

		// skip header lines
		for i := 0; i < 7; i++ {
			line, err = r.ReadString('\n')
			if err != nil {
				return err
			}
		}

		// extract data from the table
		index := 1
		for {
			line, err = r.ReadString('\n')
			if err == io.EOF {
				err = nil
				break
			}
			if err != nil {
				log.Fatalf("error parsing table [%s]: %v\n",
					fname,
					err,
				)
			}

			// handle critical or minimum ionisation energy value
			if strings.Contains(line, "M") {
				k := 0.0
				_, err = fmt.Sscanf(line, "%f", &k)
				if err != nil {
					log.Fatalf("error parsing table [%s]:\nline=%q\nerr=%v\n",
						fname,
						line,
						err,
					)
				}
				k *= 1e-3
				if strings.Contains(line, "Minimum") {
					mat.kmi = k
				} else {
					mat.kcrit = k
				}
			} else {

				var (
					k   = 0.0
					x   = 0.0
					a   = 0.0
					be  = 0.0
					tmp = 0.0
				)
				_, err = fmt.Sscanf(line,
					"%f %f %f %f %f %f %f %f %f",
					&k, &tmp, &a, &tmp, &tmp, &tmp, &be, &tmp, &x,
				)
				if err != nil {
					log.Fatalf(
						"error parsing table [%s]:\nline=%q\nerr=%v\n",
						fname,
						line,
						err,
					)
				}

				if index >= numPDG {
					return fmt.Errorf(
						"pumas: wrong format for file [%s] (#lines=%d>=%d)",
						f.Name(),
						index,
						numPDG,
					)
				}

				k *= 1e-3  // MeV        --> GeV
				a *= 1e-4  // MeV cm^2/g --> GeV m^2/kg
				be *= 1e-4 // MeV cm^2/g -->  GeV m^2/kg
				x *= 1e1   // g/cm^2     -->  kg/m^2

				// save energy loss for last entry
				if index == numPDG-1 {
					mat.amax = a
					mat.bmax = be / (k + muon.Mass)
				}

				// update table values
				{
					de := a + be
					gamma := 1.0 + k/muon.Mass
					beta := math.Sqrt(1 - 1/(gamma*gamma))

					// end-point statistics
					mat.dE[index] = de
					physK[index] = k
					mat.x[index] = x

					// weighted integrands
					mat.t[index] = 1 / (de * beta * gamma)
					mat.lp[index] = 1 / (de * math.Sqrt(k*k+2*k*muon.Mass))
				}

				index++
			}
		}
		err = f.Close()
		if err != nil {
			log.Fatalf("pumas: error closing file [%s]: %v\n", f.Name(), err)
		}

		// compute cumulative integral for the weighted time
		// Larmor's phase and scattering angle using a trapezoidal rule
		{
			var (
				time, phase float64
				dt, dp      float64
			)

			for i := 1; i < numPDG-1; i++ {
				dk := 0.5 * (physK[i] - physK[i-1])
				dt = dk * (mat.t[i] + mat.t[i-1])
				dp = dk * (mat.lp[i] + mat.lp[i-1]) * physKLarmor
				mat.t[i-1] = time
				time += dt
				mat.lp[i-1] = phase
				phase += dp
			}

			{
				hk := (physK[numPDG-1] - physK[numPDG-2]) / (physK[numPDG-2] - physK[numPDG-3])
				mat.t[numPDG-2] = time
				mat.t[numPDG-1] = time + hk*dt
				mat.lp[numPDG-2] = phase
				mat.lp[numPDG-1] = phase + hk*dp
			}

			// set the phase origin to the maximum tabulated kinetic value
			lpmax := mat.lp[numPDG-1]
			for i, v := range mat.lp {
				mat.lp[i] = lpmax - v
			}
		}

		// compute the cumulative integrals for Larmor's deflection using
		// a trapezoidal rule
		{
			var (
				x, dx [numTaylor]float64
			)
			for i := numPDG - 2; i >= 1; i-- {
				// compute deflection starting from max energy down to 0
				dk := 0.5 * (physK[i+1] - physK[i])
				p1 := mat.lp[i]
				p2 := mat.lp[i+1]
				j1 := 1.0 / mat.dE[i]
				j2 := 1.0 / mat.dE[i+1]

				f1 := 1.0
				f2 := 1.0

				for j := 0; j < numTaylor; j++ {
					mat.li[j][i+1] = x[j]
					dx[j] = dk * (f1*j1 + f2*j2)
					x[j] += dx[j]

					f1 *= p1
					f2 *= p2
				}
			}

			// extrapolate the end points
			for j := 0; j < numTaylor; j++ {
				mat.li[j][1] = x[j]
				hk := (physK[1] - physK[0]) / (physK[2] - physK[1])
				mat.li[j][0] = x[j] + hk*dx[j]
			}
		}

		// mass composition of materials, hardcoded
		{
			mat.comp = make([]componentType, mat.numElem)
			var (
				names   []string
				weights []float64
			)
			switch MaterialKind(imed) {
			case StandardRock:
				names = []string{"Rock:Standard"}
				weights = []float64{1.0}

			case DryAir:
				names = []string{"Carbon:Amorphous", "Nitrogen:Gas",
					"Oxygen:Gas", "Argon:Gas"}
				weights = []float64{0.000124, 0.755267, 0.231781, 0.012827}

			case Concrete:
				names = []string{"Hydrogen:Gas", "Carbon:Amorphous", "Oxygen:Gas", "Sodium", "Magnesium", "Aluminium", "Silicon", "Potassium", "Calcium", "Iron"}
				weights = []float64{0.010000, 0.001000, 0.529107, 0.016000, 0.002000, 0.033872, 0.337021, 0.013000, 0.044000, 0.014000}

			case Iron:
				names = []string{"Iron"}
				weights = []float64{1.0}

			case Lead:
				names = []string{"Lead"}
				weights = []float64{1.0}

			default:
				return fmt.Errorf("pumas: invalid material kind (%d)", imed)
			}
			err = mapElements(mat, names, weights)
			if err != nil {
				log.Fatalf("error mapping elements for material %d", imed)
			}

		}
	}

	log.Printf("---------------- CSDA data loaded ----------------\n")
	log.Printf("    Muon mass = %f GeV\n", muon.Mass)
	for i := range materials {
		mat := &materials[i]
		log.Printf("o Medium %d: %s\n", i+1, mat.name)
		log.Printf("    + Density            = %.3e g/cm^3\n", mat.density*1e-3)
		log.Printf("    + Ionisation minimum = %.3e GeV\n", mat.kmi)
		log.Printf("    + Critical energy    = %.3e GeV\n", mat.kcrit)
		log.Printf("    + Composition:\n")
		for _, c := range mat.comp {
			elt := elements[c.elem]
			log.Printf(
				"        . Z = %2.0f, A = %6.3f g/mol, w = %.3e\n",
				elt.z,
				elt.a,
				c.frac,
			)
		}
	}
	return err
}

// mapElements fills the material definition for a given atomic number
func mapElements(mat *materialDef, name []string, w []float64) error {
	var err error
	for i := 0; i < mat.numElem; i++ {
		iel := getElementIndex(name[i])
		mat.comp[i].elem = iel
		mat.comp[i].frac = w[i]
	}
	return err
}

// getElementIndex finds the element index for a given element name
func getElementIndex(name string) int {
	for i, e := range elements {
		if e.name == name {
			return i
		}
	}
	return -1
}

func parseDensity(line string) (float64, error) {
	var err error
	var density float64
	tokens := strings.Split(line, "density = ")
	line = strings.TrimSpace(tokens[1])
	idx := strings.Index(line, "$")
	switch idx {
	case -1:
		density, err = strconv.ParseFloat(line, 64)
		if err != nil {
			return density, err
		}
	default:
		// FIXME(sbinet) test various inputs

		// expects 'density = $ 2.0 \times10^{- 3}$\n'
		line = strings.Trim(line, "$ ")
		_, err = fmt.Sscanf(line, "%f", &density)
		if err != nil {
			return density, err
		}
		line = strings.Split(line, "{-")[1]
		line = strings.Split(line, "}")[0]
		line = strings.TrimSpace(line)
		exp := 0.0
		_, err = fmt.Sscanf(line, "%f", &exp)
		if err != nil {
			return density, err
		}
		density *= math.Pow(10, -exp)
	}

	density *= 1000 // convert to kg/m^3
	return density, err
}

// grammage returns the total grammage (in kg/m^2) seen by a muon for the given
// material and initial kinetic energy
func grammage(mat MaterialKind, kinetic float64) float64 {
	if kinetic < physK[0] {
		return 0
	}
	if kinetic >= physK[numPDG-1] {
		// constant loss model
		material := materials[mat]
		k0 := physK[numPDG-1]
		k1 := material.amax/material.bmax - muon.Mass
		return material.x[numPDG-1] + 1/material.bmax*math.Log((kinetic-k1)/(k0-k1))
	}

	return linearInterpolate(physK[:], materials[mat].x[:], kinetic)
}

// linearInterpolate implements a linear interpolation
func linearInterpolate(xs, ys []float64, x float64) float64 {
	i1 := getIndex(xs, x)
	i2 := i1 + 1

	h := (x - xs[i1]) / (xs[i2] - xs[i1])
	return ys[i1] + h*(ys[i2]-ys[i1])
}

// getIndex returns the index (into xs) for the value x, using dichotomy
func getIndex(xs []float64, x float64) int {
	if x < xs[0] {
		return -1
	}
	imax := len(xs) - 1
	if x >= xs[imax] {
		return imax
	}

	i1, _ := refineBracket(xs, x, 0, imax)
	return i1
}

func refineBracket(xs []float64, x float64, p1, p2 int) (int, int) {
	i := (p1 + p2) / 2

	if x >= xs[i] {
		p1 = i
	} else {
		p2 = i
	}

	if p2-p1 >= 2 {
		return refineBracket(xs, x, p1, p2)
	}
	return p1, p2
}

// lengthLimit returns the length limit (if any) for the propagation with
// a locator.
func (ctx *Context) lengthLimit(fwd int, ki, xi float64, mat MaterialKind, density float64) (float64, bool) {
	length := -1.0
	switch {
	case fwd != 0 && ctx.Kinetic.Min > 0:
		if ki <= ctx.Kinetic.Min {
			return length, false
		}
		length = (xi - grammage(mat, ctx.Kinetic.Min)) / density

	case fwd == 0 && ctx.Kinetic.Max > 0:
		if ki >= ctx.Kinetic.Max {
			return length, false
		}
		length = (grammage(mat, ctx.Kinetic.Max) - xi) / density
	}
	return length, true
}

// stepThrough performs a unit step through a medium.
// stepThrough returns the index of the end step medium or a negative value if
// a boundary condition was reached.
func (ctx *Context) stepThrough(locator Locator, media []Medium, length float64, mediumIndex int, charge float64, state *State, scattering, fwd int) int {

	idist := state.Distance
	e := state.Kinetic + muon.Mass
	p := math.Sqrt(e*e - muon.Mass*muon.Mass)
	b := p / e
	medium := media[mediumIndex]
	material := medium.Material
	density := medium.Density

	// compute the Larmor radius and magnetic deflection direction
	rLarmor := 0.0
	var uT fmom.Vec3
	if medium.Magnetized {
		B := medium.Magnet
		uT[0] = state.Direction[1]*B[2] - state.Direction[2]*B[1]
		uT[1] = state.Direction[2]*B[0] - state.Direction[0]*B[2]
		uT[2] = state.Direction[0]*B[1] - state.Direction[1]*B[0]
		BT := math.Sqrt(uT[0]*uT[0] + uT[1]*uT[1] + uT[2]*uT[2])
		iBT := 1.0 / BT
		rLarmor = p * iBT / physKLarmor
		uT[0] *= iBT
		uT[1] *= iBT
		uT[2] *= iBT
	}

	// randomize the next hard collision
	var (
		t         float64
		screening [maxComps]float64
		lambda    [maxComps]float64
		mu0       float64
		ihard     int
	)
	if scattering != 0 {
		// compute the mean free paths and screening factors
		var (
			invlbM  float64
			invlb1M float64
			sML     float64
			sMH     float64
		)
		for i := 0; i < materials[material].numElem; i++ {
			comp := materials[material].comp[i]
			elem := elements[comp.elem]
			s := moliereScreening(p, b, elem.z)
			lb := wentzelPath(p, b, elem.z, elem.z, s) / (density * comp.frac)
			invlb := 1.0 / lb
			invlbM += invlb
			invlb1M += wentzelTransport1(s) * invlb

			sMH += s * invlb
			sML += 1.0 / (s * lb)
			screening[i] = s
			lambda[i] = lb
		}

		// set the hard scattering mean free path
		lbM := 1.0 / invlbM
		lbH := 0.0
		{
			lbH = ctx.TransportRatio / invlb1M
			if lbH > ctx.Step.Max {
				lbH = ctx.Step.Max
			}
			lambdaK := ctx.EnergyLossRatio * state.Kinetic / (energyLoss(material, state.Kinetic) * density)
			if lbH > lambdaK {
				lbH = lambdaK
			}
			if rLarmor > 0 {
				lambdaB := rLarmor * ctx.MagneticRatio
				if lbH > lambdaB {
					lbH = lambdaB
				}
			}
		}

		// compute the hard scattering cut-off angle
		if lbM > lbH {
			lbH = lbM
			mu0 = 0
		} else {
			sM := 0.0 // effective screening factor
			if lbH > 2*lbM {
				sM = sMH * lbM // asymptotic value for lbH >> lbM
			} else {
				sM = 1 / (sML * lbM)
			}
			mu0 = sM * (lbH - lbM) / (sM*lbH + lbM) // Asymptotic value for lbH ~= lbM
			// the former asymptotic approximation yields better than 1 percent
			// relative accuracy on mu.
			// we could further improve the result with a few newtonian
			// iterations
		}
		t = -lbH * math.Log(ctx.Rand.Float64())

		// randomize the hard scatterer element
		{
			zeta := ctx.Rand.Float64()
			invlb := 0.0
			for i := 0; i < materials[material].numElem; i++ {
				ihard = i
				invlb += 1.0 / lambda[i]
				if zeta*invlbM <= invlb {
					break
				}
			}
		}

	} else {
		// disable the soft scattering
		mu0 = 0

		// set the step path
		lambdaK := ctx.EnergyLossRatio * state.Kinetic / (energyLoss(material, state.Kinetic) * density)
		t = ctx.Step.Max
		if t > lambdaK {
			t = lambdaK
		}
		if rLarmor > 0 {
			lambdaB := rLarmor * ctx.MagneticRatio
			if t > lambdaB {
				t = lambdaB
			}
		}
	}

	// do the soft scattering step
	tau := 0.0
	if mu0 == 0 {
		tau = 0.0
	} else {
		// randomize the virtual soft scattering vertex
		tau = ctx.Rand.Float64() * t

		// update the position
		{
			newIndex := ctx.updateStepPosition(locator, medium, length, mediumIndex, tau, state, fwd)
			if newIndex != mediumIndex {
				return newIndex
			}
		}

		// update the direction
		{
			ctS := 0.0
			mu1 := 0.0
			varMu := 0.0
			lb1 := 0.0
			lb2 := 0.0
			for i := 0; i < materials[material].numElem; i++ {
				si := screening[i]
				lbi := lambda[i]
				C := (si * (1 + si)) / lbi
				r := mu0 / (si + mu0)
				L := math.Log((si + mu0) / si)
				lb1 += 2 * C * (L - r)
				lb2 += 6 * C * ((1+2*si)*L - (1+2*si+mu0)*r)
			}

			if (lb1 <= 0) || (lb2 <= 0) {
				mu1 = 0
				varMu = 0
			} else {
				lb1 = 1 / lb1
				lb2 = 1 / lb2
				mu1 = 0.5 * (1 - math.Exp(-t/lb1))
				varMu = 0.25 * ((1+2*math.Exp(-t/lb2))/3 - math.Exp(-2*t/lb1))
				if varMu < 0 {
					varMu = 0
				}
			}

			// randomize the soft scattering cosine
			if mu1 > 0 {
				zeta := ctx.Rand.Float64()
				omu := 1 - mu1
				alpha := varMu / (mu1 * omu)
				alpha = (3*alpha + math.Sqrt(alpha*(alpha+8))) / (4 * (1 - alpha))
				if zeta < omu {
					ctS = 1 - 2*mu1*math.Pow(zeta/omu, alpha)
				} else {
					ctS = 2*omu*math.Pow((1-zeta)/mu1, alpha) - 1
				}
			}

			// update the direction
			if ctS > 0 {
				ctx.rotateStepDirection(ctS, state)
			}
		}
	}

	// propagate to the step end vertex
	{
		newIndex := ctx.updateStepPosition(locator, medium, length, mediumIndex, t-tau, state, fwd)

		// apply the magnetic rotation
		if rLarmor != 0 {
			u := state.Direction
			theta := 2 * float64(fwd-1) * charge * (state.Distance - idist) / rLarmor
			sin := math.Sin(theta)
			cos := math.Cos(theta)
			state.Direction[0] = cos*u[0] + sin*uT[0]
			state.Direction[1] = cos*u[1] + sin*uT[1]
			state.Direction[2] = cos*u[2] + sin*uT[2]
		}

		if newIndex != mediumIndex {
			return newIndex
		}
	}

	// apply the hard scattering
	if scattering != 0 {
		zeta := ctx.Rand.Float64()
		s := screening[ihard]
		ctH := 1 - 2*(mu0+(s+mu0)*zeta*(1-mu0)/(s+1-zeta*(1-mu0)))
		ctx.rotateStepDirection(ctH, state)
	}

	return mediumIndex
}

// updateStepPosition updates the step point position and direction.
// updateStepPosition returns the index of the end step medium or a negative
// value if a boundary condition was reached.
func (ctx *Context) updateStepPosition(locator Locator, medium Medium, length float64, mediumIndex int, step float64, state *State, fwd int) int {

	endMedium := mediumIndex
	if length > 0 {
		// check the total length boundary condition
		d := length - state.Distance
		if d <= step {
			step = d
			endMedium = -1
		}
	}
	step *= 2 * float64(fwd-1)

	if locator != nil {
		// check for a change of medium
		state.Position[0] += step * state.Direction[0]
		state.Position[1] += step * state.Direction[1]
		state.Position[2] += step * state.Direction[2]

		index := locator(ctx, state.Position)
		if index != mediumIndex {
			// locate the medium change by dichotomy
			s1 := 0.0
			s2 := -step
			for math.Abs(s1-s2) > ctx.Step.Min {
				s3 := 0.5 * (s1 + s2)
				x := state.Position[0] + s3*state.Direction[0]
				y := state.Position[1] + s3*state.Direction[1]
				z := state.Position[2] + s3*state.Direction[2]
				if locator(ctx, fmom.Vec3{x, y, z}) != index {
					s1 = s3
				} else {
					s2 = s3
				}
			}
			state.Position[0] += s1 * state.Direction[0]
			state.Position[1] += s1 * state.Direction[1]
			state.Position[2] += s1 * state.Direction[2]
			step += s1
			if endMedium >= 0 {
				endMedium = index
			}
		}
	} else {
		state.Position[0] += step * state.Direction[0]
		state.Position[1] += step * state.Direction[1]
		state.Position[2] += step * state.Direction[2]
	}

	state.Distance += math.Abs(step)

	return mediumIndex
}

// rotateStepDirection rotates the step direction using an arbitrary phase
// origin for phi.
func (ctx *Context) rotateStepDirection(ct float64, state *State) {
	// check the numerical sine
	stsq := 1 - ct*ct
	if stsq <= 0 {
		return
	}

	st := math.Sqrt(stsq)

	// select the co-vectors for the local basis
	var (
		u0 fmom.Vec3
		a  = fmom.Vec3{
			math.Abs(state.Direction[0]),
			math.Abs(state.Direction[1]),
			math.Abs(state.Direction[2]),
		}
	)
	if a[0] > a[1] {
		if a[0] > a[2] {
			nrm := 1 / math.Sqrt(state.Direction[0]*state.Direction[0]+state.Direction[2]*state.Direction[2])
			u0[0] = -nrm * state.Direction[2]
			u0[2] = +nrm * state.Direction[0]
		} else {
			nrm := 1 / math.Sqrt(state.Direction[1]*state.Direction[1]+state.Direction[2]*state.Direction[2])
			u0[1] = +nrm * state.Direction[2]
			u0[2] = -nrm * state.Direction[0]
		}
	} else {
		if a[1] > a[2] {
			nrm := 1 / math.Sqrt(state.Direction[0]*state.Direction[0]+state.Direction[1]*state.Direction[1])
			u0[0] = +nrm * state.Direction[1]
			u0[1] = -nrm * state.Direction[0]
		} else {
			nrm := 1 / math.Sqrt(state.Direction[1]*state.Direction[1]+state.Direction[2]*state.Direction[2])
			u0[1] = +nrm * state.Direction[2]
			u0[2] = -nrm * state.Direction[1]
		}
	}

	u1 := fmom.Vec3{
		u0[0]*state.Direction[2] - u0[2]*state.Direction[1],
		u0[2]*state.Direction[0] - u0[0]*state.Direction[2],
		u0[1]*state.Direction[1] - u0[1]*state.Direction[0],
	}

	// apply the rotation
	phi := 2 * math.Pi * ctx.Rand.Float64()
	cp := math.Cos(phi)
	sp := math.Sin(phi)

	state.Direction[0] = ct*state.Direction[0] + st*(cp*u0[0]+sp*u1[0])
	state.Direction[1] = ct*state.Direction[1] + st*(cp*u0[1]+sp*u1[1])
	state.Direction[2] = ct*state.Direction[2] + st*(cp*u0[2]+sp*u1[2])

	return
}

// kinetic returns the initial kinetic energy (in GeV) for the given total
// grammage and material
func kinetic(material MaterialKind, grammage float64) float64 {
	mat := materials[material]
	if grammage < mat.x[0] {
		return 0
	}

	imax := len(mat.x) - 1
	if grammage >= mat.x[imax] {
		// constant loss model
		a := mat.amax
		b := mat.bmax
		k0 := physK[imax]
		k1 := a/b - muon.Mass
		return k1 + (k0-k1)*math.Exp(b*(grammage-mat.x[imax]))
	}

	return linearInterpolate(mat.x[:], physK[:], grammage)
}

// energyLoss computes the energy loss per unit weight in GeV/(kg*m^-2)
func energyLoss(material MaterialKind, kinetic float64) float64 {
	if kinetic < physK[0] {
		return 0
	}

	mat := materials[material]
	imax := len(physK) - 1
	if kinetic >= physK[imax] {
		// constant loss model
		return mat.amax + mat.bmax*(kinetic+muon.Mass)
	}

	return linearInterpolate(physK[:], mat.dE[:], kinetic)
}

// weightedTime returns the normalized proper time (in kg/m^2) for the given
// material and kinetic energy
func weightedTime(material MaterialKind, kinetic float64) float64 {
	if kinetic < physK[0] {
		return 0
	}

	mat := materials[material]
	imax := len(physK) - 1
	if kinetic >= physK[imax] {
		// constant loss model
		a := mat.amax
		b := mat.bmax
		e0 := physK[imax] + muon.Mass
		e1 := kinetic + muon.Mass
		return mat.t[imax] + muon.Mass/a*math.Log((e1/e0)*(a+b*e0)/(a+b*e1))
	}

	return linearInterpolate(physK[:], mat.t[:], kinetic)
}

// moliereScreening returns the Moliere fit of the angular screening for the
// scattering DCS
func moliereScreening(p, beta, z float64) float64 {
	const third = 1.0 / 3.0
	d := 2.10674530e-06 * math.Pow(z, third) / p
	a := 7.29735257e-03 * z / beta
	return d * d * (1.13 + 3.76*a*a)
}

// wentzelPath returns the scattering mean free path for a Wentzel DCS
func wentzelPath(p, beta, z, a, screening float64) float64 {
	d := p * beta / z
	return a * 2.54910918e8 * screening * (1. + screening) * (d * d)
}

// wentzelTransport1 returns the 1st transport coefficient for a Wentzel DCS.
func wentzelTransport1(screening float64) float64 {
	a1 := 1.0 + screening
	return 2.0 * screening * (a1*math.Log(a1/screening) - 1.)
}

// wentzelTransport2 returns the 2nd tranport coefficient for a Wentzel DCS
func wentzelTransport2(screening float64) float64 {
	a1 := 1. + screening
	return 6. * screening * a1 * ((1.+2.*screening)*math.Log(a1/screening) - 2.)
}

func init() {
	dir := gas.MustAbs("github.com/sbinet/tmvl/pumas/data/physics/pdg/2013")
	err := loadTables(dir)
	if err != nil {
		log.Fatalf("pumas: error loading data tables: %v\n", err)
	}
}
