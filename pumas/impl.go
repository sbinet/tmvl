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
)

const (
	numPDG       = 146
	numMaterials = 5
	numElements  = 14
	numTaylor    = 8
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

func init() {
	dir := gas.MustAbs("github.com/sbinet/tmvl/pumas/data/physics/pdg/2013")
	err := loadTables(dir)
	if err != nil {
		log.Fatalf("pumas: error loading data tables: %v\n", err)
	}
}
