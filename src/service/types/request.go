package types

type Request struct {
	Wavelet   string
	Levels    uint
	IsInverse bool
	Values    []float32
}
