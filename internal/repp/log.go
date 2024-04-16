package repp

import (
	"os"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"
)

var (
	// logLevel is a configurable log level
	verboseLogging bool

	logLevel = zap.LevelEnablerFunc(func(level zapcore.Level) bool {

		// true: log message at this level
		// false: skip message at this level
		if verboseLogging {
			return level >= zapcore.DebugLevel
		} else {
			return level >= zapcore.InfoLevel
		}
	})

	// https://pkg.go.dev/go.uber.org/zap?utm_source=godoc#AtomicLevel
	l = zap.New(
		zapcore.NewCore(
			zapcore.NewConsoleEncoder(zap.NewDevelopmentEncoderConfig()),
			zapcore.Lock(os.Stderr),
			logLevel,
		),
	)

	// rlog is the default sugared logger
	rlog = l.Sugar()
)

func SetVerboseLogging() {
	verboseLogging = true
}

func isVerboseLogging() bool {
	return verboseLogging
}
