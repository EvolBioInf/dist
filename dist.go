// Package dist implements data structures and functions for reading, writing, and manipulating distance matrices in PHYLIP format.
package dist

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"math"
	"strconv"
	"strings"
	"text/tabwriter"
)

// A DistMat holds a square matrix of distances and a slice of taxon names.
type DistMat struct {
	Matrix [][]float64
	Names  []string
}

// A DistMat is read using a Scanner.
type Scanner struct {
	r   *bufio.Reader
	mat *DistMat
}

// String returns a distance matrix as a table.
func (d *DistMat) String() string {
	mat := ""
	var buf []byte
	buffer := bytes.NewBuffer(buf)
	w := new(tabwriter.Writer)
	w.Init(buffer, 1, 0, 2, ' ', 0)
	n := len(d.Matrix)
	fmt.Fprintf(w, "%d\n", n)
	w.Flush()
	for i := 0; i < n; i++ {
		fmt.Fprintf(w, "%s", d.Names[i])
		for j := 0; j < n; j++ {
			fmt.Fprintf(w, "\t%.3g", d.Matrix[i][j])
		}
		fmt.Fprint(w, "\n")
	}
	w.Flush()
	mat = buffer.String()
	return mat
}

// The method MakeSymmetrical forms the arithmetic average between corresponding entries in a DistMat.
func (d *DistMat) MakeSymmetrical() {
	n := len(d.Matrix)
	for i := 0; i < n-1; i++ {
		for j := i + 1; j < n; j++ {
			m1 := d.Matrix[i][j]
			m2 := d.Matrix[j][i]
			a := (m1 + m2) / 2.0
			d.Matrix[i][j] = a
			d.Matrix[j][i] = a
		}
	}
}

// The method Delete removes one taxon from DistMat.
func (d *DistMat) Delete(taxon int) {
	j := 0
	for i, name := range d.Names {
		if i != taxon {
			d.Names[j] = name
			j++
		}
	}
	d.Names = d.Names[:j]
	n := len(d.Matrix)
	ii := 0
	for i := 0; i < n; i++ {
		if i == taxon {
			continue
		}
		jj := 0
		for j := 0; j < n; j++ {
			if j == taxon {
				continue
			}
			d.Matrix[ii][jj] = d.Matrix[i][j]
			jj++
		}
		ii++
	}
	n--
	for i := 0; i < n; i++ {
		d.Matrix[i] = d.Matrix[i][:n]
	}
	d.Matrix = d.Matrix[:n]
}

// Method DeletePair deletes a pair of taxa.
func (d *DistMat) DeletePair(t1, t2 int) {
	j := 0
	for i, n := range d.Names {
		if i != t1 && i != t2 {
			d.Names[j] = n
			j++
		}
	}
	d.Names = d.Names[:j]
	r := 0
	for i, row := range d.Matrix {
		if i == t1 || i == t2 {
			continue
		}
		d.Matrix[r] = row
		r++
	}
	d.Matrix = d.Matrix[:r]
	for i, row := range d.Matrix {
		c := 0
		for j, col := range row {
			if j == t1 || j == t2 {
				continue
			}
			row[c] = col
			c++
		}
		d.Matrix[i] = row[:c]
	}
}

// Method Append appends a taxon consisting of a name and a slice of distances to the DistMat.
func (d *DistMat) Append(name string, data []float64) {
	n := len(d.Matrix)
	if n != len(data) {
		i := "matrix length (%d) != data length (%d)"
		log.Fatalf(i, n, len(data))
	}
	d.Names = append(d.Names, name)
	data = append(data, 0.0)
	d.Matrix = append(d.Matrix, data)
	for i := 0; i < n; i++ {
		d.Matrix[i] = append(d.Matrix[i], data[i])
	}
}

// Method Min returns the minimum entry of the distance matrix and its indexes.
func (d *DistMat) Min() (min float64, minI, minJ int) {
	min = math.MaxFloat64
	n := len(d.Names)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j {
				if d.Matrix[i][j] < min {
					min = d.Matrix[i][j]
					minI = i
					minJ = j
				}
			}
		}
	}
	return min, minI, minJ
}

// Method Max returns the maximum entry of the distance matrix and its indexes.
func (d *DistMat) Max() (min float64, maxI, maxJ int) {
	max := -math.MaxFloat64
	n := len(d.Names)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j {
				if max < d.Matrix[i][j] {
					max = d.Matrix[i][j]
					maxI = i
					maxJ = j
				}
			}
		}
	}
	return max, maxI, maxJ
}

// Scan reads input matrix by matrix.
func (s *Scanner) Scan() bool {
	var err error
	const num = "1234567890"
	l, err := s.r.ReadString('\n')
	for err == nil && strings.IndexAny(l, num) < 0 {
		l, err = s.r.ReadString('\n')
	}
	if err != nil {
		return false
	}
	l = strings.TrimRight(l, "\r\n")
	n, err := strconv.Atoi(l)
	if err != nil {
		log.Fatalf("can't convert %q", l)
	}
	s.mat = NewDistMat(n)
	for i := 0; i < n; i++ {
		l, err = s.r.ReadString('\n')
		if err != nil {
			return false
		}
		fields := strings.Fields(l)
		s.mat.Names[i] = fields[0]
		for j := 1; j <= n; j++ {
			s.mat.Matrix[i][j-1], err =
				strconv.ParseFloat(fields[j], 64)
			if err != nil {
				log.Fatalf("can't read %q", fields[j])
			}
		}

	}
	return true
}

// The method DistanceMatrix returns the last DistMat scanned.
func (s *Scanner) DistanceMatrix() *DistMat {
	return s.mat
}

// Function NewSequence returns a new n x n DistMat.
func NewDistMat(n int) *DistMat {
	d := new(DistMat)
	d.Matrix = make([][]float64, n)
	for i := 0; i < n; i++ {
		d.Matrix[i] = make([]float64, n)
	}
	d.Names = make([]string, n)
	return d
}

// The function NewScanner returns a new Scanner.
func NewScanner(r io.Reader) *Scanner {
	sc := new(Scanner)
	sc.r = bufio.NewReader(r)
	return sc
}
