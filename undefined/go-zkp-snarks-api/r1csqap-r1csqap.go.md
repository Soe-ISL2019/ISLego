---
description: f1csqap/r1csqap.go는 r1cs 단계에서 생성된 값들을 qap 값으로 변환시키는 코드
---

# r1csqap/r1csqap.go

## source code

{% tabs %}
{% tab title="r1csqap.go" %}
```go
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
```go
package r1csqap

import (
	"bytes"
	"math/big" //큰수 계산(또는 부동 소수 점수, 유리수의 정밀도 산수)을 위한 패키지

	"github.com/arnaucube/go-snark/fields"
)
```
{% endtab %}

{% tab title="struct" %}
```go
// 다항식 필드는 다항식 연산이 수행되는 유한 위의 다항식
type PolynomialField struct {
	F fields.Fq //go-snark/fields/Fq 구조체1
}
```
{% endtab %}

{% tab title="func" %}
```go
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

```go
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

```go
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

```go
// NewPolynomialField 함수는 주어진 유한체로 새로운 다항식 필드를 생성
func NewPolynomialField(f fields.Fq) PolynomialField {
	return PolynomialField{
		f,
	}
}
```

```go
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

```go
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

```go
// max 함수는 정수형 a, b 인자를 비교하여 큰 값을 반환
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
```

```go
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

```go
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

```go
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

```go
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

```go
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

```go
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

```go
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

```go
// DivisorPolynomial 함수는 2 개의 다항식이 주어지면 나눗셈(제수) 다항식을 반환
func (pf PolynomialField) DivisorPolynomial(px, z []*big.Int) []*big.Int {
	quo, _ := pf.Div(px, z)
	return quo
}
```
{% endtab %}
{% endtabs %}

## PolynomialField  

> fields.Fq 구조체는 유한체\(mod Q상의\)를 구현한 구조체 패키지이다. 때문에 해당 패키지의 모든 연산은 mod Q 상에서 이루어진다. [https://github.com/arnaucube/go-snark/blob/master/fields/fq.go](https://github.com/arnaucube/go-snark/blob/master/fields/fq.go)

```go
// 다항식 필드는 다항식 연산이 수행되는 유한 위의 다항식
type PolynomialField struct {
	F fields.Fq //go-snark/fields/Fq 구조체1
}
```

### PolynomialField 구조체 내 메소드

* Mul\(a, b \[\]_big.Int\) \[\]_big.Int 
* Div\(a, b \[\]_big.Int\) \(\[\]_big.Int, \[\]_big.Int\)_ 
* _Add\(a, b \[\]_big.Int\) \[\]_big.Int_ 
* _Sub\(a, b \[\]_big.Int\) \[\]_big.Int_ 
* _Eval\(v \[\]_big.Int, x _big.Int\)_ big.Int 
* NewPolZeroAt\(pointPos, totalPoints int, height _big.Int\) \[\]_big.Int
*  LagrangeInterpolation\(v \[\]_big.Int\) \[\]_big.Int 
* R1CSToQAP\(a, b, c \[\]\[\]_big.Int\) \(\[\]\[\]_big.Int, \[\]\[\]_big.Int, \[\]\[\]_big.Int, \[\]_big.Int\)_ 
* _CombinePolynomials\(r \[\]_big.Int, ap, bp, cp \[\]\[\]_big.Int\) \(\[\]_big.Int, \[\]_big.Int, \[\]_big.Int, \[\]_big.Int\)_
* _DivisorPolynomial\(px, z \[\]_big.Int\) \[\]\*big.Int

