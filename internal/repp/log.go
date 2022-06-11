package repp

import (
	"os"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"
)

var (
	// LogLevel is a configurable log level
	LogLevel = zap.NewAtomicLevelAt(zap.InfoLevel)

	// https://pkg.go.dev/go.uber.org/zap?utm_source=godoc#AtomicLevel
	l = zap.New(
		zapcore.NewCore(
			zapcore.NewConsoleEncoder(zap.NewDevelopmentEncoderConfig()),
			zapcore.Lock(os.Stderr),
			LogLevel,
		),
	)

	// rlog is the default sugared logger
	rlog = l.Sugar()
)
