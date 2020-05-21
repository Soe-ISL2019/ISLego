---
description: R1CS 다항식을 생성하고 witness를 계산
---

# circuitcompiler/circuit.go

방정식을 연산이 용이한 형태인 R1CS 다항식\(A~C항과 GATE로 표시\)으로 변환과정은 아래 링크 참조

{% embed url="https://medium.com/decipher-media/zero-knowledge-proof-chapter-2-deep-dive-into-zk-snarks-f8b16e1b7b4c" %}

{% tabs %}
{% tab title="circuit.go" %}
```go
package circuitcompiler

import (
	"errors"
	"math/big" //큰수 계산(또는 부동 소수 점수, 유리수의 정밀도 산수)을 위한 패키지
	"strconv"

	"github.com/arnaucube/go-snark/r1csqap" // r1csqap import
)

// Circuit is the data structure of the compiled circuit
type Circuit struct {
	NVars         int
	NPublic       int
	NSignals      int
	PrivateInputs []string
	PublicInputs  []string
	Signals       []string     // 도대체 뭘까...
	Witness       []*big.Int   // big.int형 포인터가 저장된 배열
	Constraints   []Constraint // 다항식 구조체들의 배열
	R1CS          struct {     //R1CS 다항식은 A * B = C로 이루어져있다.
		A [][]*big.Int // A항
		B [][]*big.Int // B항
		C [][]*big.Int // C항
	}
}

// Constraint is the data structure of a flat code operation
type Constraint struct { // 다항식 GATE의 구조체
	// v1 op v2 = out
	Op      string
	V1      string
	V2      string
	Out     string
	Literal string

	PrivateInputs []string // in func declaration case
	PublicInputs  []string // in func declaration case
}

// 문자열 arr에서 문자 e를 찾아 index값을 반환하는 함수
func indexInArray(arr []string, e string) int {
	for i, a := range arr {
		if a == e {
			return i
		}
	}
	return -1
}

// 숫자로 된 문자열 값의 valid 여부와 숫자로 변환한 정수값을 반환하는 함수
func isValue(a string) (bool, int) {
	v, err := strconv.Atoi(a)
	if err != nil {
		return false, 0
	}
	return true, v
}

// variable을 insert하는 역할로 보임
func insertVar(arr []*big.Int, signals []string, v string, used map[string]bool) ([]*big.Int, map[string]bool) {
	isVal, value := isValue(v)
	valueBigInt := big.NewInt(int64(value)) // value값을 큰숫자형 int로 변환
	if isVal {                              // v가 값이 있다면
		arr[0] = new(big.Int).Add(arr[0], valueBigInt) // Add()함수의 역할을 모르겠음
	} else {
		if !used[v] { // value가 값이 사용되지 않았다면
			panic(errors.New("using variable before it's set")) // 경고창 출력
		} // value 값이 사용되었다면
		arr[indexInArray(signals, v)] = new(big.Int).Add(arr[indexInArray(signals, v)], big.NewInt(int64(1))) // Add()함수의 역할을 모르겠음
	}
	return arr, used
}

// variable을 insert하는 과정중에 negative한 느낌의 역할...?
func insertVarNeg(arr []*big.Int, signals []string, v string, used map[string]bool) ([]*big.Int, map[string]bool) {
	isVal, value := isValue(v)              // v를 숫자인 value로 변환
	valueBigInt := big.NewInt(int64(value)) // value값을 큰숫자형 int로 변환
	if isVal {                              // v가 값이 있다면
		arr[0] = new(big.Int).Add(arr[0], valueBigInt) // Add()함수의 역할을 모르겠음
	} else {
		if !used[v] { // value가 값이 사용되지 않았다면
			panic(errors.New("using variable before it's set")) // 경고창 출력
		} // value 값이 사용되었다면
		arr[indexInArray(signals, v)] = new(big.Int).Add(arr[indexInArray(signals, v)], big.NewInt(int64(-1))) // Add()함수의 역할을 모르겠음
	}
	return arr, used
}

// GenerateR1CS generates the R1CS polynomials from the Circuit
// Circuit 구조체에 연결된 메서 - R1CS 생성함수
func (circ *Circuit) GenerateR1CS() ([][]*big.Int, [][]*big.Int, [][]*big.Int) {
	// from flat code to R1CS

	// big.int를 가리키는 포인터가 저장되는 2차원 배열
	var a [][]*big.Int
	var b [][]*big.Int
	var c [][]*big.Int

	used := make(map[string]bool)                 // string-bool 형식의 해시테이블 생성, used는 string이 사용되었는지 확인하는 데이터를 가짐
	for _, constraint := range circ.Constraints { // circuit 구조체에서 Constraints 배열을 갖고와서 index는 버리고 value만 사용
		aConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals)) // ricsqap의 ArrayOfBigZeros함수 호출
		bConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals))
		cConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals))

		// if existInArray(constraint.Out) {
		// if used[constraint.Out] {
		// panic(errors.New("out variable already used: " + constraint.Out))
		// }
		used[constraint.Out] = true // 다항식 GATE 구조체의 결과값 out이 사용되었다고 표시
		if constraint.Op == "in" {  // 다항식 GATE 구조체의 연산자가 "in"이면
			for i := 0; i <= len(circ.PublicInputs); i++ {
				aConstraint[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Add(aConstraint[indexInArray(circ.Signals, constraint.Out)], big.NewInt(int64(1)))
				aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.Out, used)
				bConstraint[0] = big.NewInt(int64(1))
			}
			continue

		} else if constraint.Op == "+" { // 다항식 GATE 구조체의 연산자가 + 이면
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V1, used)
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V2, used)
			bConstraint[0] = big.NewInt(int64(1))
		} else if constraint.Op == "-" { // 다항식 GATE 구조체의 연산자가 - 이면
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVarNeg(aConstraint, circ.Signals, constraint.V1, used)
			aConstraint, used = insertVarNeg(aConstraint, circ.Signals, constraint.V2, used)
			bConstraint[0] = big.NewInt(int64(1))
		} else if constraint.Op == "*" { // 다항식 GATE 구조체의 연산자가 * 이면
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V1, used)
			bConstraint, used = insertVar(bConstraint, circ.Signals, constraint.V2, used)
		} else if constraint.Op == "/" { // 다항식 GATE 구조체의 연산자가 / 이면
			cConstraint, used = insertVar(cConstraint, circ.Signals, constraint.V1, used)
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			bConstraint, used = insertVar(bConstraint, circ.Signals, constraint.V2, used)
		}

		a = append(a, aConstraint) // 2차원 배열 a에 aConstraint를 추가
		b = append(b, bConstraint)
		c = append(c, cConstraint)

	}
	// 다항식 GATE의 구조체에 2차원 배열 a를 저장
	circ.R1CS.A = a
	circ.R1CS.B = b
	circ.R1CS.C = c
	return a, b, c
}

func grabVar(signals []string, w []*big.Int, vStr string) *big.Int {
	isVal, v := isValue(vStr)    // vStr 값의 valid 여부를 반환하고 값을 정수로 변환
	vBig := big.NewInt(int64(v)) // 값 -> int형 -> big.int형으로 변환
	if isVal {                   // vStr값이 valid하면
		return vBig // 변환한 값을 return
	} else { // vStr값이 invalid하면
		return w[indexInArray(signals, vStr)] // signals 문자열에서 vStr값에 대한 index를 포인터배열 w의 index로 하여 값을 return
	}
}

// input 구조체?
type Inputs struct {
	Private []*big.Int
	Public  []*big.Int
}

// CalculateWitness calculates the Witness of a Circuit based on the given inputs
// witness = [ one, output, publicInputs, privateInputs, ...]
func (circ *Circuit) CalculateWitness(privateInputs []*big.Int, publicInputs []*big.Int) ([]*big.Int, error) {
	if len(privateInputs) != len(circ.PrivateInputs) {
		return []*big.Int{}, errors.New("given privateInputs != circuit.PublicInputs")
	}
	if len(publicInputs) != len(circ.PublicInputs) {
		return []*big.Int{}, errors.New("given publicInputs != circuit.PublicInputs")
	}
	w := r1csqap.ArrayOfBigZeros(len(circ.Signals)) // ricsqap의 ArrayOfBigZeros함수 호출
	w[0] = big.NewInt(int64(1))                     // 1 -> int형 -> big.Int형으로 변환하고 w배열의 0번째 index에 저장
	for i, input := range publicInputs {
		w[i+1] = input
	}
	for i, input := range privateInputs {
		w[i+len(publicInputs)+1] = input
	}
	for _, constraint := range circ.Constraints {
		if constraint.Op == "in" {
		} else if constraint.Op == "+" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Add(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "-" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Sub(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "*" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Mul(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "/" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Div(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		}
	}
	return w, nil
}

```
{% endtab %}

{% tab title="import" %}
```go
import (
	"errors"
	"math/big" //큰수 계산(또는 부동 소수 점수, 유리수의 정밀도 산수)을 위한 패키지
	"strconv"

	"github.com/arnaucube/go-snark/r1csqap" // r1csqap import
)
```
{% endtab %}

{% tab title="struct" %}
```go
type Circuit struct {
	NVars         int
	NPublic       int
	NSignals      int
	PrivateInputs []string
	PublicInputs  []string
	Signals       []string     // 도대체 뭘까...
	Witness       []*big.Int   // big.int형 포인터가 저장된 배열
	Constraints   []Constraint // 다항식 구조체들의 배열
	R1CS          struct {     //R1CS 다항식은 A * B = C로 이루어져있다.
		A [][]*big.Int // A항
		B [][]*big.Int // B항
		C [][]*big.Int // C항
	}
}

// Constraint is the data structure of a flat code operation
type Constraint struct { // 다항식 GATE의 구조체
	// v1 op v2 = out
	Op      string
	V1      string
	V2      string
	Out     string
	Literal string

	PrivateInputs []string // in func declaration case
	PublicInputs  []string // in func declaration case
}
	
// input 구조체?
type Inputs struct {
	Private []*big.Int
	Public  []*big.Int
}
```
{% endtab %}

{% tab title="function" %}
```go
// 문자열 arr에서 문자 e를 찾아 index값을 반환하는 함수
func indexInArray(arr []string, e string) int {
	for i, a := range arr {
		if a == e {
			return i
		}
	}
	return -1
}

// 숫자로 된 문자열 값의 valid 여부와 숫자로 변환한 정수값을 반환하는 함수
func isValue(a string) (bool, int) {
	v, err := strconv.Atoi(a)
	if err != nil {
		return false, 0
	}
	return true, v
}

// variable을 insert하는 역할로 보임
func insertVar(arr []*big.Int, signals []string, v string, used map[string]bool) ([]*big.Int, map[string]bool) {
	isVal, value := isValue(v)
	valueBigInt := big.NewInt(int64(value)) // value값을 큰숫자형 int로 변환
	if isVal {                              // v가 값이 있다면
		arr[0] = new(big.Int).Add(arr[0], valueBigInt) // Add()함수의 역할을 모르겠음
	} else {
		if !used[v] { // value가 값이 사용되지 않았다면
			panic(errors.New("using variable before it's set")) // 경고창 출력
		} // value 값이 사용되었다면
		arr[indexInArray(signals, v)] = new(big.Int).Add(arr[indexInArray(signals, v)], big.NewInt(int64(1))) // Add()함수의 역할을 모르겠음
	}
	return arr, used
}

// variable을 insert하는 과정중에 negative한 느낌의 역할...?
func insertVarNeg(arr []*big.Int, signals []string, v string, used map[string]bool) ([]*big.Int, map[string]bool) {
	isVal, value := isValue(v)              // v를 숫자인 value로 변환
	valueBigInt := big.NewInt(int64(value)) // value값을 큰숫자형 int로 변환
	if isVal {                              // v가 값이 있다면
		arr[0] = new(big.Int).Add(arr[0], valueBigInt) // Add()함수의 역할을 모르겠음
	} else {
		if !used[v] { // value가 값이 사용되지 않았다면
			panic(errors.New("using variable before it's set")) // 경고창 출력
		} // value 값이 사용되었다면
		arr[indexInArray(signals, v)] = new(big.Int).Add(arr[indexInArray(signals, v)], big.NewInt(int64(-1))) // Add()함수의 역할을 모르겠음
	}
	return arr, used
}

func grabVar(signals []string, w []*big.Int, vStr string) *big.Int {
	isVal, v := isValue(vStr)    // vStr 값의 valid 여부를 반환하고 값을 정수로 변환
	vBig := big.NewInt(int64(v)) // 값 -> int형 -> big.int형으로 변환
	if isVal {                   // vStr값이 valid하면
		return vBig // 변환한 값을 return
	} else { // vStr값이 invalid하면
		return w[indexInArray(signals, vStr)] // signals 문자열에서 vStr값에 대한 index를 포인터배열 w의 index로 하여 값을 return
	}
}
```
{% endtab %}

{% tab title="method" %}
```go
// GenerateR1CS generates the R1CS polynomials from the Circuit
// Circuit 구조체에 연결된 메서드 - R1CS 생성함수
func (circ *Circuit) GenerateR1CS() ([][]*big.Int, [][]*big.Int, [][]*big.Int) {
	// from flat code to R1CS

	// big.int를 가리키는 포인터가 저장되는 2차원 배열
	var a [][]*big.Int
	var b [][]*big.Int
	var c [][]*big.Int

	used := make(map[string]bool)                 // string-bool 형식의 해시테이블 생성, used는 string이 사용되었는지 확인하는 데이터를 가짐
	for _, constraint := range circ.Constraints { // circuit 구조체에서 Constraints 배열을 갖고와서 index는 버리고 value만 사용
		aConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals)) // ricsqap의 ArrayOfBigZeros함수 호출
		bConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals))
		cConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals))

		// if existInArray(constraint.Out) {
		// if used[constraint.Out] {
		// panic(errors.New("out variable already used: " + constraint.Out))
		// }
		used[constraint.Out] = true // 다항식 GATE 구조체의 결과값 out이 사용되었다고 표시
		if constraint.Op == "in" {  // 다항식 GATE 구조체의 연산자가 "in"이면
			for i := 0; i <= len(circ.PublicInputs); i++ {
				aConstraint[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Add(aConstraint[indexInArray(circ.Signals, constraint.Out)], big.NewInt(int64(1)))
				aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.Out, used)
				bConstraint[0] = big.NewInt(int64(1))
			}
			continue

		} else if constraint.Op == "+" { // 다항식 GATE 구조체의 연산자가 + 이면
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V1, used)
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V2, used)
			bConstraint[0] = big.NewInt(int64(1))
		} else if constraint.Op == "-" { // 다항식 GATE 구조체의 연산자가 - 이면
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVarNeg(aConstraint, circ.Signals, constraint.V1, used)
			aConstraint, used = insertVarNeg(aConstraint, circ.Signals, constraint.V2, used)
			bConstraint[0] = big.NewInt(int64(1))
		} else if constraint.Op == "*" { // 다항식 GATE 구조체의 연산자가 * 이면
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V1, used)
			bConstraint, used = insertVar(bConstraint, circ.Signals, constraint.V2, used)
		} else if constraint.Op == "/" { // 다항식 GATE 구조체의 연산자가 / 이면
			cConstraint, used = insertVar(cConstraint, circ.Signals, constraint.V1, used)
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			bConstraint, used = insertVar(bConstraint, circ.Signals, constraint.V2, used)
		}

		a = append(a, aConstraint) // 2차원 배열 a에 aConstraint를 추가
		b = append(b, bConstraint)
		c = append(c, cConstraint)

	}
	// 다항식 GATE의 구조체에 2차원 배열 a를 저장
	circ.R1CS.A = a
	circ.R1CS.B = b
	circ.R1CS.C = c
	return a, b, c
}

// CalculateWitness calculates the Witness of a Circuit based on the given inputs
// witness = [ one, output, publicInputs, privateInputs, ...]
func (circ *Circuit) CalculateWitness(privateInputs []*big.Int, publicInputs []*big.Int) ([]*big.Int, error) {
	if len(privateInputs) != len(circ.PrivateInputs) {
		return []*big.Int{}, errors.New("given privateInputs != circuit.PublicInputs")
	}
	if len(publicInputs) != len(circ.PublicInputs) {
		return []*big.Int{}, errors.New("given publicInputs != circuit.PublicInputs")
	}
	w := r1csqap.ArrayOfBigZeros(len(circ.Signals)) // ricsqap의 ArrayOfBigZeros함수 호출
	w[0] = big.NewInt(int64(1))                     // 1 -> int형 -> big.Int형으로 변환하고 w배열의 0번째 index에 저장
	for i, input := range publicInputs {
		w[i+1] = input
	}
	for i, input := range privateInputs {
		w[i+len(publicInputs)+1] = input
	}
	for _, constraint := range circ.Constraints {
		if constraint.Op == "in" {
		} else if constraint.Op == "+" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Add(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "-" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Sub(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "*" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Mul(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "/" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Div(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		}
	}
	return w, nil
}
```
{% endtab %}
{% endtabs %}



