package repp

import (
	"fmt"
	"math"
)

type seqScores struct {
	gcContent            float64
	longestHomopolymer   int
	min50WindowGCContent float64
	max50WindowGCContent float64
}

func (s seqScores) String() string {
	return fmt.Sprintf("gc=%f, hp=%d, min50gc=%f, max50gc=%f",
		s.gcContent, s.longestHomopolymer, s.min50WindowGCContent, s.max50WindowGCContent)
}

func (s *seqScores) add(s1 seqScores) {
	s.gcContent = math.Max(s.gcContent, s1.gcContent)
	s.longestHomopolymer = int(math.Max(float64(s.longestHomopolymer), float64(s1.longestHomopolymer)))
	s.min50WindowGCContent = math.Max(s.min50WindowGCContent, s1.min50WindowGCContent)
	s.max50WindowGCContent = math.Max(s.max50WindowGCContent, s1.max50WindowGCContent)
}

type scoreAlg interface {
	nextBp(bp rune, pos int)
	score() float64
}

type gcScore struct {
	seqLen  int
	gcCount int
}

func (a *gcScore) nextBp(bp rune, pos int) {
	if bp == 'C' || bp == 'G' {
		a.gcCount++
	}
}

func (a *gcScore) score() float64 {
	return float64(a.gcCount) / float64(a.seqLen)
}

type windowGCScore struct {
	windowSeq             []rune
	windowLen             int
	windowGCCount         int
	windowGCCountBoundary int
	boundaryCmp           func(v1, v2 int) bool
}

func (a *windowGCScore) nextBp(bp rune, pos int) {
	if bp == 'C' || bp == 'G' {
		a.windowGCCount++
	}
	if pos < a.windowLen {
		a.windowSeq[pos] = bp
		a.windowGCCountBoundary = a.windowGCCount
	} else {
		prevWindowSeq := a.windowSeq
		a.windowSeq = append(prevWindowSeq[1:], bp)
		if prevWindowSeq[0] == 'C' || prevWindowSeq[0] == 'G' {
			a.windowGCCount--
		}
		if a.boundaryCmp(a.windowGCCount, a.windowGCCountBoundary) {
			a.windowGCCountBoundary = a.windowGCCount
		}
	}
}

func (a *windowGCScore) score() float64 {
	return float64(a.windowGCCountBoundary) / float64(a.windowLen)
}

type homopolymerScore struct {
	prevBp                 rune
	currentHomopolymerSize int
	longestHomopolymer     int
}

func (a *homopolymerScore) nextBp(bp rune, pos int) {
	if bp == a.prevBp {
		a.currentHomopolymerSize++
	} else {
		if a.currentHomopolymerSize > a.longestHomopolymer {
			a.longestHomopolymer = a.currentHomopolymerSize
		}
		a.prevBp = bp
		a.currentHomopolymerSize = 1
	}
}

func (a *homopolymerScore) score() float64 {
	if a.currentHomopolymerSize > a.longestHomopolymer {
		return float64(a.currentHomopolymerSize)
	}
	return float64(a.longestHomopolymer)
}

func fragSeqQualityChecks(seq string) seqScores {

	gcContent := &gcScore{
		seqLen:  len(seq),
		gcCount: 0,
	}
	homopolymerCount := &homopolymerScore{
		prevBp:                 rune(0),
		currentHomopolymerSize: 0,
		longestHomopolymer:     0,
	}
	minWindowGCContent := &windowGCScore{
		windowSeq:             make([]rune, 50),
		windowLen:             50,
		windowGCCount:         0,
		windowGCCountBoundary: 0,
		boundaryCmp:           func(v1, v2 int) bool { return v1 < v2 },
	}
	maxWindowGCContent := &windowGCScore{
		windowSeq:             make([]rune, 50),
		windowLen:             50,
		windowGCCount:         0,
		windowGCCountBoundary: 0,
		boundaryCmp:           func(v1, v2 int) bool { return v1 > v2 },
	}
	algs := []scoreAlg{
		gcContent,
		homopolymerCount,
		minWindowGCContent,
		maxWindowGCContent,
	}

	for pos, bp := range seq {
		for _, alg := range algs {
			alg.nextBp(bp, pos)
		}
	}

	return seqScores{
		gcContent:            gcContent.score(),
		longestHomopolymer:   int(homopolymerCount.score()),
		min50WindowGCContent: minWindowGCContent.score(),
		max50WindowGCContent: maxWindowGCContent.score(),
	}
}
