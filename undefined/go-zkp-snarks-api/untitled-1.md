---
description: >-
  bn128/g1.go의 G1 구조체는 타원곡선의 연산을 위한 객체로, Jacobian Coordinates에서 타원곡선의 연산을 하는
  구조체라고 생각됨
---

# bn128/g1.go 분석

## 포함 패키지 

```go
import (
	"math/big" //큰수 계산(또는 부동 소수 점수, 유리수의 정밀도 산수)을 위한 패키지
	"github.com/arnaucube/go-snark/fields" //실제 go snark git에 존재하는 패키지 => 체(Field)를 구현하고 있음
)
```

## G1 구조체 

> fields.Fq 구조체는 유한체\(mod Q상의\)를 구현한 구조체 패키지이다. 때문에 해당 패키지의 모든 연산은 mod Q 상에서 이루어진다. [https://github.com/arnaucube/go-snark/blob/master/fields/fq.go](https://github.com/arnaucube/go-snark/blob/master/fields/fq.go)
>
> big.Int 구조체는 매우 큰 정수 연산을 위해 존재하는 패키지이다. 때문에 간단한 사칙연산 부터 모든 연산들이 메소드로 구현되어 있다. [https://golang.org/pkg/math/big/\#Int](https://golang.org/pkg/math/big/#Int)
>
> 결과적으로 G1은 타원곡선 알고리즘을 구현하고 있다. [https://en.wikipedia.org/wiki/Elliptic\_curve\_point\_multiplication\#Double-and-add](https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#Double-and-add)

```go

//야코비안 프라임 곡선, y^2 = x^3 + ax + b 에서 타원 곡선 점을 나타내는데 사용
type G1 struct { //G1 구조체 역할은 단지 관련 메소드를 사용하기 위함
	F fields.Fq   //go-snark/fields/Fq 구조체1
	G [3]*big.Int // 큰 Int 포인터 배열
}
```

### G1 구조체 내 메소드 

* NewG1\(f fields.Fq, g \[2\]\*big.Int\) G1
* Zero\(\) \[2\]\*big.Int
* IsZero\(p \[3\]\*big.Int\) bool
* Add\(p1, p2 \[3\]_big.Int\) \[3\]_big.Int
* Neg\(p \[3\]big.Int\) \[3\]big.Int
* Sub\(a, b \[3\]big.Int\) \[3\]big.Int
* Double\(p \[3\]big.Int\) \[3\]big.Int
* MulScalar\(p \[3\]big.Int, e big.Int\) \[3\]\*big.Int
* Affine\(p \[3\]big.Int\) \[2\]big.Int
* Equal\(p1, p2 \[3\]\*big.Int\) bool

## G1.NewG1\(f fields.Fq, g \[2\]\*big.Int\) G1

> 매개변수로 받아온 field\(체\)구조체와 x,y 좌표를 받고, z좌표를 1로 초기화하여, G1 구조체를 생성한다.

```go
func NewG1(f fields.Fq, g [2]*big.Int) G1 {
	var g1 G1
	g1.F = f
	g1.G = [3]*big.Int{ //큰 int 포인 배열 3개 할당
		g[0], // 매개변수 배열 g
		g[1],
		g1.F.One(), //field.Fq 구조체의 One()함수 호출 => 1로 선언된 int 생성
		/* field.Fq.One() 함수
		func (fq Fq) One() *big.Int {
			return big.NewInt(int64(1)) // 1로 선언된 int 생성
		}*/
	}
	return g1
}
```

## G1.Zero\(\) \[2\]\*big.Int

> 0으로 초기화된 크기 2의 big.Int 배열을 할당해줌

```go
func (g1 G1) Zero() [2]*big.Int {
	return [2]*big.Int{g1.F.Zero(), g1.F.Zero()} //g1.F는 field.Fq 구조체, .Zero() 매소드는 .One()과 같으며, 0으로 선언된 int 생성해줌
	/* field.Fq.Zero() 함수
	func (fq Fq) Zero() *big.Int {
		return big.NewInt(int64(0)) // 1로 선언된 int 생성
	}*/
}
```

## IsZero\(p \[3\]\*big.Int\) bool

> 주어진 점의 p\[2\] 값, 즉 z 축 값이 0인지 확인

```go
func (g1 G1) IsZero(p [3]*big.Int) bool {
	return g1.F.IsZero(p[2]) //p[2]가 0인지 확인, 여기서 p[2]는 z좌표, G1구조체에서는 NewOne()으로 초기화 해줬음
}
```

## Add\(p1, p2 \[3\]big.Int\) \[3\]big.Int

> 타원 곡선 위의 두 점의 덧셈연산을 하는 메소드,

{% code title="알고리즘" %}
```text
 U1 = X1*Z2^2
 U2 = X2*Z1^2
 S1 = Y1*Z2^3
 S2 = Y2*Z1^3
 if (U1 == U2)
   if (S1 != S2)
     return POINT_AT_INFINITY
   else 
     return POINT_DOUBLE(X1, Y1, Z1)
 H = U2 - U1
 R = S2 - S1
 X3 = R^2 - H^3 - 2*U1*H^2
 Y3 = R*(U1*H^2 - X3) - S1*H^3
 Z3 = H*Z1*Z2
 return (X3, Y3, Z3)
```
{% endcode %}

```go
func (g1 G1) Add(p1, p2 [3]*big.Int) [3]*big.Int { // point Addition (12M + 4S)

	// https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates
	// https://github.com/zcash/zcash/blob/master/src/snark/libsnark/algebra/curves/alt_bn128/alt_bn128_g1.cpp#L208
	// http://hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2007-bl.op3

	if g1.IsZero(p1) { //p1[2]가 0인지 확인
		return p2
	}
	if g1.IsZero(p2) { //p2[2]가 0인지 확인
		return p1
	}

	// p1, p2 배열은 점의 좌표를 나타냄
	x1 := p1[0] //
	y1 := p1[1] //
	z1 := p1[2] //
	x2 := p2[0] //
	y2 := p2[1] //
	z2 := p2[2] //

	z1z1 := g1.F.Square(z1) //g1.F => field.Fq 구조체, Fq.Square() 메소드 호출 Z1Z1 = Z1^2
	z2z2 := g1.F.Square(z2) //말 그대로 제곱 연산 Z2Z2 = Z2^2
	/*
		func (fq Fq) Square(a *big.Int) *big.Int {
			m := new(big.Int).Mul(a, a)	//a*a연산
			return new(big.Int).Mod(m, fq.Q)
		}*/

	u1 := g1.F.Mul(x1, z2z2) //U1 = X1*Z2Z2
	u2 := g1.F.Mul(x2, z1z1) //U2 = X2*Z1Z1
	/*
		func (fq Fq) Mul(a, b *big.Int) *big.Int {
			m := new(big.Int).Mul(a, b)
			return new(big.Int).Mod(m, fq.Q)
		}*/

	t0 := g1.F.Mul(z2, z2z2) //t0 = Z2*Z2Z2			t0 = Z2^3
	s1 := g1.F.Mul(y1, t0)   //S1 = Y1*t0				S1 = Y1*Z2^3

	t1 := g1.F.Mul(z1, z1z1) //t1 = Z1*Z1Z1			t1 = Z1^3
	s2 := g1.F.Mul(y2, t1)   //S2 = Y2*t1				S2 = Y2*Z1^3

	h := g1.F.Sub(u2, u1)  //H = U2-U1					H = X2*Z1^2-X1*Z2^2
	t2 := g1.F.Add(h, h)   //t2 = H+H = H*2			t2 = 2*(X2*Z1^2-X1*Z2^2)
	i := g1.F.Square(t2)   //I = t2^2					I = 4*(X2*Z1^2-X1*Z2^2)^2
	j := g1.F.Mul(h, i)    //J = H*I					J = 4*(X2*Z1^2-X1*Z2^2)^3
	t3 := g1.F.Sub(s2, s1) //t3 = S2-S1				t3 = Y2*Z1^3-Y1*Z2^3
	r := g1.F.Add(t3, t3)  //r = 2*t3					r = 2*(Y2*Z1^3-Y1*Z2^3)
	v := g1.F.Mul(u1, i)   //V = U1*I					V = X1*Z2^2*4*(X2*Z1^2-X1*Z2^2)^2
	t4 := g1.F.Square(r)   //t4 = r^2					t4 = 4*(Y2*Z1^3-Y1*Z2^3)^2
	t5 := g1.F.Add(v, v)   //t5 = 2*V
	t6 := g1.F.Sub(t4, j)  //t6 = t4-J
	x3 := g1.F.Sub(t6, t5) //X3 = t6-t5
	t7 := g1.F.Sub(v, x3)  //t7 = V-X3
	t8 := g1.F.Mul(s1, j)  //t8 = S1*J
	t9 := g1.F.Add(t8, t8) //t9 = 2*t8
	t10 := g1.F.Mul(r, t7) //t10 = r*t7

	y3 := g1.F.Sub(t10, t9) //Y3 = t10-t9

	t11 := g1.F.Add(z1, z2)    //t11 = Z1+Z2
	t12 := g1.F.Square(t11)    //t12 = t11^2
	t13 := g1.F.Sub(t12, z1z1) //t13 = t12-Z1Z1
	t14 := g1.F.Sub(t13, z2z2) //t14 = t13-Z2Z2
	z3 := g1.F.Mul(t14, h)     //Z3 = t14*H

	return [3]*big.Int{x3, y3, z3}
}
```

## Neg\(p \[3\]big.Int\) \[3\]big.Int

> 타원 곡선 상의 덧셈의 역원이 되는 점을 반환, 즉 {p\[0\], -p\[1\], p\[2\]}을 반환하게 됨

```go
func (g1 G1) Neg(p [3]*big.Int) [3]*big.Int {
	return [3]*big.Int{
		p[0],
		g1.F.Neg(p[1]), //p[1]의 덧셈의 역원을 구함, 주로 y1으로 사용되므로, -y1
		p[2],
	}
}
```

## Sub\(a, b \[3\]big.Int\) \[3\]big.Int

> 타원 곡선에서 덧셈의 역원이 되는 점을 Add\(\)함수를 통해 더해줌. 즉 뺄셈 연산

```go
func (g1 G1) Sub(a, b [3]*big.Int) [3]*big.Int {
	return g1.Add(a, g1.Neg(b)) //즉 p[1]만 덧셈의 역원으로 하여 Add 함수 
}
```

## Double\(p \[3\]big.Int\) \[3\]big.Int

> 타원 곡선 위에서 점 P의 2배 연산을 해주는 함수

```go
func (g1 G1) Double(p [3]*big.Int) [3]*big.Int { //point Doubling (4M + 6S or 4M + 4S)

	// https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates
	// http://hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/doubling/dbl-2009-l.op3
	// https://github.com/zcash/zcash/blob/master/src/snark/libsnark/algebra/curves/alt_bn128/alt_bn128_g1.cpp#L325

	if g1.IsZero(p) {
		return p
	}

	a := g1.F.Square(p[0]) //A = X1^2
	b := g1.F.Square(p[1]) //B = Y1^2
	c := g1.F.Square(b)    //C = B^2

	t0 := g1.F.Add(p[0], b) //t0 = X1+B
	t1 := g1.F.Square(t0)   //t1 = t0^2
	t2 := g1.F.Sub(t1, a)   //t2 = t1-A
	t3 := g1.F.Sub(t2, c)   //t3 = t2-C

	d := g1.F.Double(t3)             //D = 2*t3
	e := g1.F.Add(g1.F.Add(a, a), a) //E = 3*A
	f := g1.F.Square(e)              //F = E^2

	t4 := g1.F.Double(d)  //t4 = 2*D
	x3 := g1.F.Sub(f, t4) //X3 = F-t4

	t5 := g1.F.Sub(d, x3)         //t5 = D-X3
	twoC := g1.F.Add(c, c)        //2*C
	fourC := g1.F.Add(twoC, twoC) //4*C
	t6 := g1.F.Add(fourC, fourC)  //t6 = 8*C
	t7 := g1.F.Mul(e, t5)         //t7 = E*t5
	y3 := g1.F.Sub(t7, t6)        //Y3 = t7-t6

	t8 := g1.F.Mul(p[1], p[2]) //t8 = Y1*Z1
	z3 := g1.F.Double(t8)      //Z3 = 2*t8

	return [3]*big.Int{x3, y3, z3}
}
```

## MulScalar\(p \[3\]\*bih.Int, e \*big.Int\) \[3\]\*big.Int

> 매개변수 e의 비트열을 갖고, 타원 곡선 위의 점간 곱셈연산을 해준다. 타원 곡선위 곱연산은 다음과 같은 구조를



