# \[2020.05.21.Thurs\]

## \[cli - G. J. RA\]

## \[Groth16, r1csqap - T. H. KIM\]

### r1csqap.go

{% tabs %}
{% tab title="r1csqap.go" %}
```text
package r1csqap

import (
	"bytes"
	"math/big"

	"github.com/arnaucube/go-snark/fields"
)

// Transpose transposes the *big.Int matrix
func Transpose(matrix [][]*big.Int) [][]*big.Int {
	var r [][]*big.Int
	for i := 0; i < len(matrix[0]); i++ {
		var row []*big.Int
		for j := 0; j < len(matrix); j++ {
			row = append(row, matrix[j][i])
		}
		r = append(r, row)
	}
	return r
}

// ArrayOfBigZeros creates a *big.Int array with n elements to zero
func ArrayOfBigZeros(num int) []*big.Int {
	bigZero := big.NewInt(int64(0))
	var r []*big.Int
	for i := 0; i < num; i++ {
		r = append(r, bigZero)
	}
	return r
}
func BigArraysEqual(a, b []*big.Int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !bytes.Equal(a[i].Bytes(), b[i].Bytes()) {
			return false
		}
	}
	return true
}

// PolynomialField is the Polynomial over a Finite Field where the polynomial operations are performed
type PolynomialField struct {
	F fields.Fq
}

// NewPolynomialField creates a new PolynomialField with the given FiniteField
func NewPolynomialField(f fields.Fq) PolynomialField {
	return PolynomialField{
		f,
	}
}

// Mul multiplies two polinomials over the Finite Field
func (pf PolynomialField) Mul(a, b []*big.Int) []*big.Int {
	r := ArrayOfBigZeros(len(a) + len(b) - 1)
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(b); j++ {
			r[i+j] = pf.F.Add(
				r[i+j],
				pf.F.Mul(a[i], b[j]))
		}
	}
	return r
}

// Div divides two polinomials over the Finite Field, returning the result and the remainder
func (pf PolynomialField) Div(a, b []*big.Int) ([]*big.Int, []*big.Int) {
	// https://en.wikipedia.org/wiki/Division_algorithm
	r := ArrayOfBigZeros(len(a) - len(b) + 1)
	rem := a
	for len(rem) >= len(b) {
		l := pf.F.Div(rem[len(rem)-1], b[len(b)-1])
		pos := len(rem) - len(b)
		r[pos] = l
		aux := ArrayOfBigZeros(pos)
		aux1 := append(aux, l)
		aux2 := pf.Sub(rem, pf.Mul(b, aux1))
		rem = aux2[:len(aux2)-1]
	}
	return r, rem
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// Add adds two polinomials over the Finite Field
func (pf PolynomialField) Add(a, b []*big.Int) []*big.Int {
	r := ArrayOfBigZeros(max(len(a), len(b)))
	for i := 0; i < len(a); i++ {
		r[i] = pf.F.Add(r[i], a[i])
	}
	for i := 0; i < len(b); i++ {
		r[i] = pf.F.Add(r[i], b[i])
	}
	return r
}

// Sub subtracts two polinomials over the Finite Field
func (pf PolynomialField) Sub(a, b []*big.Int) []*big.Int {
	r := ArrayOfBigZeros(max(len(a), len(b)))
	for i := 0; i < len(a); i++ {
		r[i] = pf.F.Add(r[i], a[i])
	}
	for i := 0; i < len(b); i++ {
		r[i] = pf.F.Sub(r[i], b[i])
	}
	return r
}

// Eval evaluates the polinomial over the Finite Field at the given value x
func (pf PolynomialField) Eval(v []*big.Int, x *big.Int) *big.Int {
	r := big.NewInt(int64(0))
	for i := 0; i < len(v); i++ {
		xi := pf.F.Exp(x, big.NewInt(int64(i)))
		elem := pf.F.Mul(v[i], xi)
		r = pf.F.Add(r, elem)
	}
	return r
}

// NewPolZeroAt generates a new polynomial that has value zero at the given value
func (pf PolynomialField) NewPolZeroAt(pointPos, totalPoints int, height *big.Int) []*big.Int {
	fac := 1
	for i := 1; i < totalPoints+1; i++ {
		if i != pointPos {
			fac = fac * (pointPos - i)
		}
	}
	facBig := big.NewInt(int64(fac))
	hf := pf.F.Div(height, facBig)
	r := []*big.Int{hf}
	for i := 1; i < totalPoints+1; i++ {
		if i != pointPos {
			ineg := big.NewInt(int64(-i))
			b1 := big.NewInt(int64(1))
			r = pf.Mul(r, []*big.Int{ineg, b1})
		}
	}
	return r
}

// LagrangeInterpolation performs the Lagrange Interpolation / Lagrange Polynomials operation
func (pf PolynomialField) LagrangeInterpolation(v []*big.Int) []*big.Int {
	// https://en.wikipedia.org/wiki/Lagrange_polynomial
	var r []*big.Int
	for i := 0; i < len(v); i++ {
		r = pf.Add(r, pf.NewPolZeroAt(i+1, len(v), v[i]))
	}
	//
	return r
}

// R1CSToQAP converts the R1CS values to the QAP values
func (pf PolynomialField) R1CSToQAP(a, b, c [][]*big.Int) ([][]*big.Int, [][]*big.Int, [][]*big.Int, []*big.Int) {
	aT := Transpose(a)
	bT := Transpose(b)
	cT := Transpose(c)
	var alphas [][]*big.Int
	for i := 0; i < len(aT); i++ {
		alphas = append(alphas, pf.LagrangeInterpolation(aT[i]))
	}
	var betas [][]*big.Int
	for i := 0; i < len(bT); i++ {
		betas = append(betas, pf.LagrangeInterpolation(bT[i]))
	}
	var gammas [][]*big.Int
	for i := 0; i < len(cT); i++ {
		gammas = append(gammas, pf.LagrangeInterpolation(cT[i]))
	}
	z := []*big.Int{big.NewInt(int64(1))}
	for i := 1; i < len(alphas)-1; i++ {
		z = pf.Mul(
			z,
			[]*big.Int{
				pf.F.Neg(
					big.NewInt(int64(i))),
				big.NewInt(int64(1)),
			})
	}
	return alphas, betas, gammas, z
}

// CombinePolynomials combine the given polynomials arrays into one, also returns the P(x)
func (pf PolynomialField) CombinePolynomials(r []*big.Int, ap, bp, cp [][]*big.Int) ([]*big.Int, []*big.Int, []*big.Int, []*big.Int) {
	var ax []*big.Int
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*big.Int{r[i]}, ap[i])
		ax = pf.Add(ax, m)
	}
	var bx []*big.Int
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*big.Int{r[i]}, bp[i])
		bx = pf.Add(bx, m)
	}
	var cx []*big.Int
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*big.Int{r[i]}, cp[i])
		cx = pf.Add(cx, m)
	}

	px := pf.Sub(pf.Mul(ax, bx), cx)
	return ax, bx, cx, px
}

// DivisorPolynomial returns the divisor polynomial given two polynomials
func (pf PolynomialField) DivisorPolynomial(px, z []*big.Int) []*big.Int {
	quo, _ := pf.Div(px, z)
	return quo
}

```
{% endtab %}

{% tab title="import" %}
```
package r1csqap

import (
	"bytes"
	"math/big"

	"github.com/arnaucube/go-snark/fields"
)
```
{% endtab %}

{% tab title="struct" %}
```
// 다항식 필드는 다항식 연산이 수행되는 유한 위의 다항식
type PolynomialField struct {
	F fields.Fq
}
```
{% endtab %}

{% tab title="func" %}
```
// Transpose 함수는 *big.Int 행렬을 전치
func Transpose(matrix [][]*big.Int) [][]*big.Int {
	var r [][]*big.Int
	for i := 0; i < len(matrix[0]); i++ {
		var row []*big.Int
		for j := 0; j < len(matrix); j++ {
			row = append(row, matrix[j][i])
		}
		r = append(r, row)
	}
	return r
}
```

```text
// ArrayOfBigZeros 함는 n 개의 요소가 0 인 *big.Int 배열을 생성
func ArrayOfBigZeros(num int) []*big.Int {
	bigZero := big.NewInt(int64(0))
	var r []*big.Int
	for i := 0; i < num; i++ {
		r = append(r, bigZero)
	}
	return r
}
```

```text
// BigArraysEqual 함수는 2 개의 *big.Int 배열 길이 비교 후 같을 경우 값이 일치한지 비
func BigArraysEqual(a, b []*big.Int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !bytes.Equal(a[i].Bytes(), b[i].Bytes()) {
			return false
		}
	}
	return true
}
```

```text
// NewPolynomialField 함수는 주어진 유한체로 새로운 다항식 필드를 생성
func NewPolynomialField(f fields.Fq) PolynomialField {
	return PolynomialField{
		f,
	}
}
```

```text
// Mul 함수는 유한체에서 두 개의 다항식을 곱셈 연산
func (pf PolynomialField) Mul(a, b []*big.Int) []*big.Int {
	r := ArrayOfBigZeros(len(a) + len(b) - 1)
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(b); j++ {
			r[i+j] = pf.F.Add(
				r[i+j],
				pf.F.Mul(a[i], b[j]))
		}
	}
	return r
}
```

```text
// Div 함수는 유한체에서 두 개의 다항식을 나누고 결과와 나머지를 반환
func (pf PolynomialField) Div(a, b []*big.Int) ([]*big.Int, []*big.Int) {
	// https://en.wikipedia.org/wiki/Division_algorithm
	r := ArrayOfBigZeros(len(a) - len(b) + 1)
	rem := a
	for len(rem) >= len(b) {
		l := pf.F.Div(rem[len(rem)-1], b[len(b)-1])
		pos := len(rem) - len(b)
		r[pos] = l
		aux := ArrayOfBigZeros(pos)
		aux1 := append(aux, l)
		aux2 := pf.Sub(rem, pf.Mul(b, aux1))
		rem = aux2[:len(aux2)-1]
	}
	return r, rem
}
```

```text
// max 함수는 정수형 a, b 인자를 비교하여 큰 값을 반환
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
```

```text
// Add 함수는 유한체 위에 두 개의 다항식을 덧셈 연산
func (pf PolynomialField) Add(a, b []*big.Int) []*big.Int {
	r := ArrayOfBigZeros(max(len(a), len(b)))
	for i := 0; i < len(a); i++ {
		r[i] = pf.F.Add(r[i], a[i])
	}
	for i := 0; i < len(b); i++ {
		r[i] = pf.F.Add(r[i], b[i])
	}
	return r
}
```

```text
// Sub 유한체에서 두 개의 다항식을 뺄셈 연산 
func (pf PolynomialField) Sub(a, b []*big.Int) []*big.Int {
	r := ArrayOfBigZeros(max(len(a), len(b)))
	for i := 0; i < len(a); i++ {
		r[i] = pf.F.Add(r[i], a[i])
	}
	for i := 0; i < len(b); i++ {
		r[i] = pf.F.Sub(r[i], b[i])
	}
	return r
}
```

```text
// Eval 주어진 값 x에서 유한체에서의 다항식을 평가
func (pf PolynomialField) Eval(v []*big.Int, x *big.Int) *big.Int {
	r := big.NewInt(int64(0))
	for i := 0; i < len(v); i++ {
		xi := pf.F.Exp(x, big.NewInt(int64(i)))
		elem := pf.F.Mul(v[i], xi)
		r = pf.F.Add(r, elem)
	}
	return r
}
```

```text
// NewPolZeroAt 함수는 주어진 값에서 값이 0 인 새로운 다항식을 생성
func (pf PolynomialField) NewPolZeroAt(pointPos, totalPoints int, height *big.Int) []*big.Int {
	fac := 1
	for i := 1; i < totalPoints+1; i++ {
		if i != pointPos {
			fac = fac * (pointPos - i)
		}
	}
	facBig := big.NewInt(int64(fac))
	hf := pf.F.Div(height, facBig)
	r := []*big.Int{hf}
	for i := 1; i < totalPoints+1; i++ {
		if i != pointPos {
			ineg := big.NewInt(int64(-i))
			b1 := big.NewInt(int64(1))
			r = pf.Mul(r, []*big.Int{ineg, b1})
		}
	}
	return r
}
```

```text
// LagrangeInterpolation 함수는 라그랑주 보간법/ 라그랑주 다항식 연산을 수행
func (pf PolynomialField) LagrangeInterpolation(v []*big.Int) []*big.Int {
	// https://en.wikipedia.org/wiki/Lagrange_polynomial
	var r []*big.Int
	for i := 0; i < len(v); i++ {
		r = pf.Add(r, pf.NewPolZeroAt(i+1, len(v), v[i]))
	}
	//
	return r
}
```

```text
// R1CSToQAP 함수는 R1CS 값을 QAP 값으로 변환
func (pf PolynomialField) R1CSToQAP(a, b, c [][]*big.Int) ([][]*big.Int, [][]*big.Int, [][]*big.Int, []*big.Int) {
	aT := Transpose(a)
	bT := Transpose(b)
	cT := Transpose(c)
	var alphas [][]*big.Int
	for i := 0; i < len(aT); i++ {
		alphas = append(alphas, pf.LagrangeInterpolation(aT[i]))
	}
	var betas [][]*big.Int
	for i := 0; i < len(bT); i++ {
		betas = append(betas, pf.LagrangeInterpolation(bT[i]))
	}
	var gammas [][]*big.Int
	for i := 0; i < len(cT); i++ {
		gammas = append(gammas, pf.LagrangeInterpolation(cT[i]))
	}
	z := []*big.Int{big.NewInt(int64(1))}
	for i := 1; i < len(alphas)-1; i++ {
		z = pf.Mul(
			z,
			[]*big.Int{
				pf.F.Neg(
					big.NewInt(int64(i))),
				big.NewInt(int64(1)),
			})
	}
	return alphas, betas, gammas, z
}
```

```text
// CombinePolynomials 함수는 주어진 다항식 배열을 하나로 결합하고, P(x)를 반환
func (pf PolynomialField) CombinePolynomials(r []*big.Int, ap, bp, cp [][]*big.Int) ([]*big.Int, []*big.Int, []*big.Int, []*big.Int) {
	var ax []*big.Int
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*big.Int{r[i]}, ap[i])
		ax = pf.Add(ax, m)
	}
	var bx []*big.Int
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*big.Int{r[i]}, bp[i])
		bx = pf.Add(bx, m)
	}
	var cx []*big.Int
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*big.Int{r[i]}, cp[i])
		cx = pf.Add(cx, m)
	}

	px := pf.Sub(pf.Mul(ax, bx), cx)
	return ax, bx, cx, px
}
```

```text
// DivisorPolynomial 함수는 2 개의 다항식이 주어지면 나눗셈(제수) 다항식을 반환
func (pf PolynomialField) DivisorPolynomial(px, z []*big.Int) []*big.Int {
	quo, _ := pf.Div(px, z)
	return quo
}
```
{% endtab %}
{% endtabs %}

### groth16.go

{% tabs %}
{% tab title="groth16.go" %}
```
// implementation of https://eprint.iacr.org/2016/260.pdf

package groth16

import (
	"fmt"
	"math/big"

	"github.com/arnaucube/go-snark/bn128"
	"github.com/arnaucube/go-snark/circuitcompiler"
	"github.com/arnaucube/go-snark/fields"
	"github.com/arnaucube/go-snark/r1csqap"
)

type Pk struct { // Proving Key
	BACDelta [][3]*big.Int // {( βui(x)+αvi(x)+wi(x) ) / δ } from l+1 to m
	Z        []*big.Int
	G1       struct {
		Alpha    [3]*big.Int
		Beta     [3]*big.Int
		Delta    [3]*big.Int
		At       [][3]*big.Int // {a(τ)} from 0 to m
		BACGamma [][3]*big.Int // {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m
	}
	G2 struct {
		Beta     [3][2]*big.Int
		Gamma    [3][2]*big.Int
		Delta    [3][2]*big.Int
		BACGamma [][3][2]*big.Int // {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m
	}
	PowersTauDelta [][3]*big.Int // powers of τ encrypted in G1 curve, divided by δ
}
type Vk struct {
	IC [][3]*big.Int
	G1 struct {
		Alpha [3]*big.Int
	}
	G2 struct {
		Beta  [3][2]*big.Int
		Gamma [3][2]*big.Int
		Delta [3][2]*big.Int
	}
}

// Setup is the data structure holding the Trusted Setup data. The Setup.Toxic sub struct must be destroyed after the GenerateTrustedSetup function is completed
type Setup struct {
	Toxic struct {
		T      *big.Int // trusted setup secret
		Kalpha *big.Int
		Kbeta  *big.Int
		Kgamma *big.Int
		Kdelta *big.Int
	}

	// public
	Pk Pk
	Vk Vk
}

// Proof contains the parameters to proof the zkSNARK
type Proof struct {
	PiA [3]*big.Int
	PiB [3][2]*big.Int
	PiC [3]*big.Int
}

type utils struct {
	Bn  bn128.Bn128
	FqR fields.Fq
	PF  r1csqap.PolynomialField
}

// Utils is the data structure holding the BN128, FqR Finite Field over R, PolynomialField, that will be used inside the snarks operations
var Utils = prepareUtils()

func prepareUtils() utils {
	bn, err := bn128.NewBn128()
	if err != nil {
		panic(err)
	}
	// new Finite Field
	fqR := fields.NewFq(bn.R)
	// new Polynomial Field
	pf := r1csqap.NewPolynomialField(fqR)

	return utils{
		Bn:  bn,
		FqR: fqR,
		PF:  pf,
	}
}

// GenerateTrustedSetup generates the Trusted Setup from a compiled Circuit. The Setup.Toxic sub data structure must be destroyed
func GenerateTrustedSetup(witnessLength int, circuit circuitcompiler.Circuit, alphas, betas, gammas [][]*big.Int) (Setup, error) {
	var setup Setup
	var err error

	// generate random t value
	setup.Toxic.T, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}

	setup.Toxic.Kalpha, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kbeta, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kgamma, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kdelta, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}

	// z pol
	zpol := []*big.Int{big.NewInt(int64(1))}
	for i := 1; i < len(alphas)-1; i++ {
		zpol = Utils.PF.Mul(
			zpol,
			[]*big.Int{
				Utils.FqR.Neg(
					big.NewInt(int64(i))),
				big.NewInt(int64(1)),
			})
	}
	setup.Pk.Z = zpol
	zt := Utils.PF.Eval(zpol, setup.Toxic.T)
	invDelta := Utils.FqR.Inverse(setup.Toxic.Kdelta)
	ztinvDelta := Utils.FqR.Mul(invDelta, zt)

	// encrypt t values with curve generators
	// powers of tau divided by delta
	var ptd [][3]*big.Int
	ini := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, ztinvDelta)
	ptd = append(ptd, ini)
	tEncr := setup.Toxic.T
	for i := 1; i < len(zpol); i++ {
		ptd = append(ptd, Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, Utils.FqR.Mul(tEncr, ztinvDelta)))
		tEncr = Utils.FqR.Mul(tEncr, setup.Toxic.T)
	}
	// powers of τ encrypted in G1 curve, divided by δ
	// (G1 * τ) / δ
	setup.Pk.PowersTauDelta = ptd

	setup.Pk.G1.Alpha = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kalpha)
	setup.Pk.G1.Beta = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kbeta)
	setup.Pk.G1.Delta = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kdelta)
	setup.Pk.G2.Beta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kbeta)
	setup.Pk.G2.Delta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kdelta)

	setup.Vk.G1.Alpha = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kalpha)
	setup.Vk.G2.Beta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kbeta)
	setup.Vk.G2.Gamma = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kgamma)
	setup.Vk.G2.Delta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kdelta)

	for i := 0; i < len(circuit.Signals); i++ {
		// Pk.G1.At: {a(τ)} from 0 to m
		at := Utils.PF.Eval(alphas[i], setup.Toxic.T)
		a := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, at)
		setup.Pk.G1.At = append(setup.Pk.G1.At, a)

		bt := Utils.PF.Eval(betas[i], setup.Toxic.T)
		g1bt := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, bt)
		g2bt := Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, bt)
		// G1.BACGamma: {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m in G1
		setup.Pk.G1.BACGamma = append(setup.Pk.G1.BACGamma, g1bt)
		// G2.BACGamma: {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m in G2
		setup.Pk.G2.BACGamma = append(setup.Pk.G2.BACGamma, g2bt)
	}

	zero3 := [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	for i := 0; i < circuit.NPublic+1; i++ {
		setup.Pk.BACDelta = append(setup.Pk.BACDelta, zero3)
	}
	for i := circuit.NPublic + 1; i < circuit.NVars; i++ {
		// TODO calculate all at, bt, ct outside, to avoid repeating calculations
		at := Utils.PF.Eval(alphas[i], setup.Toxic.T)
		bt := Utils.PF.Eval(betas[i], setup.Toxic.T)
		ct := Utils.PF.Eval(gammas[i], setup.Toxic.T)
		c := Utils.FqR.Mul(
			invDelta,
			Utils.FqR.Add(
				Utils.FqR.Add(
					Utils.FqR.Mul(at, setup.Toxic.Kbeta),
					Utils.FqR.Mul(bt, setup.Toxic.Kalpha),
				),
				ct,
			),
		)
		g1c := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, c)

		// Pk.BACDelta: {( βui(x)+αvi(x)+wi(x) ) / δ } from l+1 to m
		setup.Pk.BACDelta = append(setup.Pk.BACDelta, g1c)
	}

	for i := 0; i <= circuit.NPublic; i++ {
		at := Utils.PF.Eval(alphas[i], setup.Toxic.T)
		bt := Utils.PF.Eval(betas[i], setup.Toxic.T)
		ct := Utils.PF.Eval(gammas[i], setup.Toxic.T)
		ic := Utils.FqR.Mul(
			Utils.FqR.Inverse(setup.Toxic.Kgamma),
			Utils.FqR.Add(
				Utils.FqR.Add(
					Utils.FqR.Mul(at, setup.Toxic.Kbeta),
					Utils.FqR.Mul(bt, setup.Toxic.Kalpha),
				),
				ct,
			),
		)
		g1ic := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, ic)
		// used in verifier
		setup.Vk.IC = append(setup.Vk.IC, g1ic)
	}

	return setup, nil
}

// GenerateProofs generates all the parameters to proof the zkSNARK from the Circuit, Setup and the Witness
func GenerateProofs(circuit circuitcompiler.Circuit, pk Pk, w []*big.Int, px []*big.Int) (Proof, error) {
	var proof Proof
	proof.PiA = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiB = Utils.Bn.Fq6.Zero()
	proof.PiC = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}

	r, err := Utils.FqR.Rand()
	if err != nil {
		return Proof{}, err
	}
	s, err := Utils.FqR.Rand()
	if err != nil {
		return Proof{}, err
	}

	// piBG1 will hold all the same than proof.PiB but in G1 curve
	piBG1 := [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}

	for i := 0; i < circuit.NVars; i++ {
		proof.PiA = Utils.Bn.G1.Add(proof.PiA, Utils.Bn.G1.MulScalar(pk.G1.At[i], w[i]))
		piBG1 = Utils.Bn.G1.Add(piBG1, Utils.Bn.G1.MulScalar(pk.G1.BACGamma[i], w[i]))
		proof.PiB = Utils.Bn.G2.Add(proof.PiB, Utils.Bn.G2.MulScalar(pk.G2.BACGamma[i], w[i]))
	}
	for i := circuit.NPublic + 1; i < circuit.NVars; i++ {
		proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(pk.BACDelta[i], w[i]))
	}

	// piA = (Σ from 0 to m (pk.A * w[i])) + pk.Alpha1 + r * δ
	proof.PiA = Utils.Bn.G1.Add(proof.PiA, pk.G1.Alpha)
	deltaR := Utils.Bn.G1.MulScalar(pk.G1.Delta, r)
	proof.PiA = Utils.Bn.G1.Add(proof.PiA, deltaR)

	// piBG1 = (Σ from 0 to m (pk.B1 * w[i])) + pk.g1.Beta + s * δ
	// piB = piB2 = (Σ from 0 to m (pk.B2 * w[i])) + pk.g2.Beta + s * δ
	piBG1 = Utils.Bn.G1.Add(piBG1, pk.G1.Beta)
	proof.PiB = Utils.Bn.G2.Add(proof.PiB, pk.G2.Beta)
	deltaSG1 := Utils.Bn.G1.MulScalar(pk.G1.Delta, s)
	piBG1 = Utils.Bn.G1.Add(piBG1, deltaSG1)
	deltaSG2 := Utils.Bn.G2.MulScalar(pk.G2.Delta, s)
	proof.PiB = Utils.Bn.G2.Add(proof.PiB, deltaSG2)

	hx := Utils.PF.DivisorPolynomial(px, pk.Z) // maybe move this calculation to a previous step

	// piC = (Σ from l+1 to m (w[i] * (pk.g1.Beta + pk.g1.Alpha + pk.C)) + h(tau)) / δ) + piA*s + r*piB - r*s*δ
	for i := 0; i < len(hx); i++ {
		proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(pk.PowersTauDelta[i], hx[i]))
	}
	proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(proof.PiA, s))
	proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(piBG1, r))
	negRS := Utils.FqR.Neg(Utils.FqR.Mul(r, s))
	proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(pk.G1.Delta, negRS))

	return proof, nil
}

// VerifyProof verifies over the BN128 the Pairings of the Proof
func VerifyProof(vk Vk, proof Proof, publicSignals []*big.Int, debug bool) bool {

	icPubl := vk.IC[0]
	for i := 0; i < len(publicSignals); i++ {
		icPubl = Utils.Bn.G1.Add(icPubl, Utils.Bn.G1.MulScalar(vk.IC[i+1], publicSignals[i]))
	}

	if !Utils.Bn.Fq12.Equal(
		Utils.Bn.Pairing(proof.PiA, proof.PiB),
		Utils.Bn.Fq12.Mul(
			Utils.Bn.Pairing(vk.G1.Alpha, vk.G2.Beta),
			Utils.Bn.Fq12.Mul(
				Utils.Bn.Pairing(icPubl, vk.G2.Gamma),
				Utils.Bn.Pairing(proof.PiC, vk.G2.Delta)))) {
		if debug {
			fmt.Println("❌ groth16 verification not passed")
		}
		return false
	}
	if debug {
		fmt.Println("✓ groth16 verification passed")
	}

	return true
}

```
{% endtab %}

{% tab title="import" %}
```go
// implementation of https://eprint.iacr.org/2016/260.pdf

package groth16

import (
	"fmt"
	"math/big"

	"github.com/arnaucube/go-snark/bn128"
	"github.com/arnaucube/go-snark/circuitcompiler"
	"github.com/arnaucube/go-snark/fields"
	"github.com/arnaucube/go-snark/r1csqap"
)
```
{% endtab %}

{% tab title="Struct" %}
```
type Pk struct { // Proving Key
	BACDelta [][3]*big.Int //  l+1 {( βui(x)+αvi(x)+wi(x) ) / δ } 
	Z        []*big.Int
	G1       struct {
		Alpha    [3]*big.Int
		Beta     [3]*big.Int
		Delta    [3]*big.Int
		At       [][3]*big.Int // {a(τ)} from 0 to m
		BACGamma [][3]*big.Int // {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m
	}
	G2 struct {
		Beta     [3][2]*big.Int
		Gamma    [3][2]*big.Int
		Delta    [3][2]*big.Int
		BACGamma [][3][2]*big.Int // {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m
	}
	PowersTauDelta [][3]*big.Int // powers of τ encrypted in G1 curve, divided by δ
}
type Vk struct { // Verification Key
	IC [][3]*big.Int
	G1 struct {
		Alpha [3]*big.Int
	}
	G2 struct {
		Beta  [3][2]*big.Int
		Gamma [3][2]*big.Int
		Delta [3][2]*big.Int
	}
}
```

```text
// Setup 구조체는 신뢰할 수있는 설정 데이터를 보유하는 데이터 구조 
// GenerateTrustedSetup 함수가 완료된 후 Setup.Toxic 하위 구조체를 삭제해야 함 
type Setup struct {
	Toxic struct {
		T      *big.Int // trusted setup secret
		Kalpha *big.Int
		Kbeta  *big.Int
		Kgamma *big.Int
		Kdelta *big.Int
	}

	// public
	Pk Pk
	Vk Vk
}

// Proof 구조체는 zk-SNARK를 증명하는 매개 변수를 포함
type Proof struct {
	PiA [3]*big.Int
	PiB [3][2]*big.Int
	PiC [3]*big.Int
}

type utils struct {
	Bn  bn128.Bn128
	FqR fields.Fq
	PF  r1csqap.PolynomialFie
```
{% endtab %}

{% tab title="func" %}
```
// Utils is the data structure holding the BN128, FqR Finite Field over R, PolynomialField, that will be used inside the snarks operations
var Utils = prepareUtils()

func prepareUtils() utils {
	bn, err := bn128.NewBn128()
	if err != nil {
		panic(err)
	}
	// new Finite Field
	fqR := fields.NewFq(bn.R)
	// new Polynomial Field
	pf := r1csqap.NewPolynomialField(fqR)

	return utils{
		Bn:  bn,
		FqR: fqR,
		PF:  pf,
	}
}
```

```text
// GenerateTrustedSetup generates the Trusted Setup from a compiled Circuit. The Setup.Toxic sub data structure must be destroyed
func GenerateTrustedSetup(witnessLength int, circuit circuitcompiler.Circuit, alphas, betas, gammas [][]*big.Int) (Setup, error) {
	var setup Setup
	var err error

	// generate random t value
	setup.Toxic.T, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}

	setup.Toxic.Kalpha, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kbeta, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kgamma, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kdelta, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}

	// z pol
	zpol := []*big.Int{big.NewInt(int64(1))}
	for i := 1; i < len(alphas)-1; i++ {
		zpol = Utils.PF.Mul(
			zpol,
			[]*big.Int{
				Utils.FqR.Neg(
					big.NewInt(int64(i))),
				big.NewInt(int64(1)),
			})
	}
	setup.Pk.Z = zpol
	zt := Utils.PF.Eval(zpol, setup.Toxic.T)
	invDelta := Utils.FqR.Inverse(setup.Toxic.Kdelta)
	ztinvDelta := Utils.FqR.Mul(invDelta, zt)

	// encrypt t values with curve generators
	// powers of tau divided by delta
	var ptd [][3]*big.Int
	ini := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, ztinvDelta)
	ptd = append(ptd, ini)
	tEncr := setup.Toxic.T
	for i := 1; i < len(zpol); i++ {
		ptd = append(ptd, Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, Utils.FqR.Mul(tEncr, ztinvDelta)))
		tEncr = Utils.FqR.Mul(tEncr, setup.Toxic.T)
	}
	// powers of τ encrypted in G1 curve, divided by δ
	// (G1 * τ) / δ
	setup.Pk.PowersTauDelta = ptd

	setup.Pk.G1.Alpha = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kalpha)
	setup.Pk.G1.Beta = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kbeta)
	setup.Pk.G1.Delta = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kdelta)
	setup.Pk.G2.Beta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kbeta)
	setup.Pk.G2.Delta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kdelta)

	setup.Vk.G1.Alpha = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kalpha)
	setup.Vk.G2.Beta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kbeta)
	setup.Vk.G2.Gamma = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kgamma)
	setup.Vk.G2.Delta = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kdelta)

	for i := 0; i < len(circuit.Signals); i++ {
		// Pk.G1.At: {a(τ)} from 0 to m
		at := Utils.PF.Eval(alphas[i], setup.Toxic.T)
		a := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, at)
		setup.Pk.G1.At = append(setup.Pk.G1.At, a)

		bt := Utils.PF.Eval(betas[i], setup.Toxic.T)
		g1bt := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, bt)
		g2bt := Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, bt)
		// G1.BACGamma: {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m in G1
		setup.Pk.G1.BACGamma = append(setup.Pk.G1.BACGamma, g1bt)
		// G2.BACGamma: {( βui(x)+αvi(x)+wi(x) ) / γ } from 0 to m in G2
		setup.Pk.G2.BACGamma = append(setup.Pk.G2.BACGamma, g2bt)
	}

	zero3 := [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	for i := 0; i < circuit.NPublic+1; i++ {
		setup.Pk.BACDelta = append(setup.Pk.BACDelta, zero3)
	}
	for i := circuit.NPublic + 1; i < circuit.NVars; i++ {
		// TODO calculate all at, bt, ct outside, to avoid repeating calculations
		at := Utils.PF.Eval(alphas[i], setup.Toxic.T)
		bt := Utils.PF.Eval(betas[i], setup.Toxic.T)
		ct := Utils.PF.Eval(gammas[i], setup.Toxic.T)
		c := Utils.FqR.Mul(
			invDelta,
			Utils.FqR.Add(
				Utils.FqR.Add(
					Utils.FqR.Mul(at, setup.Toxic.Kbeta),
					Utils.FqR.Mul(bt, setup.Toxic.Kalpha),
				),
				ct,
			),
		)
		g1c := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, c)

		// Pk.BACDelta: {( βui(x)+αvi(x)+wi(x) ) / δ } from l+1 to m
		setup.Pk.BACDelta = append(setup.Pk.BACDelta, g1c)
	}

	for i := 0; i <= circuit.NPublic; i++ {
		at := Utils.PF.Eval(alphas[i], setup.Toxic.T)
		bt := Utils.PF.Eval(betas[i], setup.Toxic.T)
		ct := Utils.PF.Eval(gammas[i], setup.Toxic.T)
		ic := Utils.FqR.Mul(
			Utils.FqR.Inverse(setup.Toxic.Kgamma),
			Utils.FqR.Add(
				Utils.FqR.Add(
					Utils.FqR.Mul(at, setup.Toxic.Kbeta),
					Utils.FqR.Mul(bt, setup.Toxic.Kalpha),
				),
				ct,
			),
		)
		g1ic := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, ic)
		// used in verifier
		setup.Vk.IC = append(setup.Vk.IC, g1ic)
	}

	return setup, nil
}
```

```text
// GenerateProofs generates all the parameters to proof the zkSNARK from the Circuit, Setup and the Witness
func GenerateProofs(circuit circuitcompiler.Circuit, pk Pk, w []*big.Int, px []*big.Int) (Proof, error) {
	var proof Proof
	proof.PiA = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiB = Utils.Bn.Fq6.Zero()
	proof.PiC = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}

	r, err := Utils.FqR.Rand()
	if err != nil {
		return Proof{}, err
	}
	s, err := Utils.FqR.Rand()
	if err != nil {
		return Proof{}, err
	}

	// piBG1 will hold all the same than proof.PiB but in G1 curve
	piBG1 := [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}

	for i := 0; i < circuit.NVars; i++ {
		proof.PiA = Utils.Bn.G1.Add(proof.PiA, Utils.Bn.G1.MulScalar(pk.G1.At[i], w[i]))
		piBG1 = Utils.Bn.G1.Add(piBG1, Utils.Bn.G1.MulScalar(pk.G1.BACGamma[i], w[i]))
		proof.PiB = Utils.Bn.G2.Add(proof.PiB, Utils.Bn.G2.MulScalar(pk.G2.BACGamma[i], w[i]))
	}
	for i := circuit.NPublic + 1; i < circuit.NVars; i++ {
		proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(pk.BACDelta[i], w[i]))
	}

	// piA = (Σ from 0 to m (pk.A * w[i])) + pk.Alpha1 + r * δ
	proof.PiA = Utils.Bn.G1.Add(proof.PiA, pk.G1.Alpha)
	deltaR := Utils.Bn.G1.MulScalar(pk.G1.Delta, r)
	proof.PiA = Utils.Bn.G1.Add(proof.PiA, deltaR)

	// piBG1 = (Σ from 0 to m (pk.B1 * w[i])) + pk.g1.Beta + s * δ
	// piB = piB2 = (Σ from 0 to m (pk.B2 * w[i])) + pk.g2.Beta + s * δ
	piBG1 = Utils.Bn.G1.Add(piBG1, pk.G1.Beta)
	proof.PiB = Utils.Bn.G2.Add(proof.PiB, pk.G2.Beta)
	deltaSG1 := Utils.Bn.G1.MulScalar(pk.G1.Delta, s)
	piBG1 = Utils.Bn.G1.Add(piBG1, deltaSG1)
	deltaSG2 := Utils.Bn.G2.MulScalar(pk.G2.Delta, s)
	proof.PiB = Utils.Bn.G2.Add(proof.PiB, deltaSG2)

	hx := Utils.PF.DivisorPolynomial(px, pk.Z) // maybe move this calculation to a previous step

	// piC = (Σ from l+1 to m (w[i] * (pk.g1.Beta + pk.g1.Alpha + pk.C)) + h(tau)) / δ) + piA*s + r*piB - r*s*δ
	for i := 0; i < len(hx); i++ {
		proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(pk.PowersTauDelta[i], hx[i]))
	}
	proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(proof.PiA, s))
	proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(piBG1, r))
	negRS := Utils.FqR.Neg(Utils.FqR.Mul(r, s))
	proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(pk.G1.Delta, negRS))

	return proof, nil
}
```

```text
// VerifyProof verifies over the BN128 the Pairings of the Proof
func VerifyProof(vk Vk, proof Proof, publicSignals []*big.Int, debug bool) bool {

	icPubl := vk.IC[0]
	for i := 0; i < len(publicSignals); i++ {
		icPubl = Utils.Bn.G1.Add(icPubl, Utils.Bn.G1.MulScalar(vk.IC[i+1], publicSignals[i]))
	}

	if !Utils.Bn.Fq12.Equal(
		Utils.Bn.Pairing(proof.PiA, proof.PiB),
		Utils.Bn.Fq12.Mul(
			Utils.Bn.Pairing(vk.G1.Alpha, vk.G2.Beta),
			Utils.Bn.Fq12.Mul(
				Utils.Bn.Pairing(icPubl, vk.G2.Gamma),
				Utils.Bn.Pairing(proof.PiC, vk.G2.Delta)))) {
		if debug {
			fmt.Println("❌ groth16 verification not passed")
		}
		return false
	}
	if debug {
		fmt.Println("✓ groth16 verification passed")
	}

	return true
}
```
{% endtab %}
{% endtabs %}

## \[bn128, circuitcompiler - D. H. GWAK, H. J. SONG, N. J. LEE\]

### bn128/g1.go

{% page-ref page="untitled-1.md" %}



