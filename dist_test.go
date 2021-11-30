package dist

import (
	"io/ioutil"
	"os"
	"testing"
)

func checkMat(d *DistMat, file string, t *testing.T) {
	res, err := ioutil.ReadFile(file)
	if err != nil {
		t.Errorf("can't open %q", file)
	}
	want := string(res)
	get := d.String()
	if get != want {
		t.Errorf("get:\n%s\nwant:\n%s", get, want)
	}
}
func TestDist(t *testing.T) {
	fn := "data/dm.phy"
	f, err := os.Open(fn)
	if err != nil {
		t.Errorf("can't open %q", fn)
	}
	defer f.Close()
	sc := NewScanner(f)
	sc.Scan()
	dm := sc.DistanceMatrix()
	checkMat(dm, "data/r1.txt", t)
	dm.MakeSymmetrical()
	checkMat(dm, "data/r2.txt", t)
	dm.Delete(0)
	checkMat(dm, "data/r3.txt", t)
	data := []float64{0.1, 0.1, 0.1, 0.1}
	dm.Append("new", data)
	checkMat(dm, "data/r4.txt", t)
	m, i, j := dm.Min()
	if m != 0.0489 || i != 0 || j != 3 {
		t.Errorf("can't find minimum")
	}
	dm.DeletePair(0, 4)
	checkMat(dm, "data/r5.txt", t)
}
