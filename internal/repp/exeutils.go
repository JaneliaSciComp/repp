package repp

import (
	"os"
)

func getExecutable(exeHomeEnvVar, binSubDir, exeName string) string {
	exeHome := os.Getenv(exeHomeEnvVar)
	if exeHome == "" {
		return exeName
	}
	if binSubDir == "" {
		return exeHome + "/" + exeName
	} else {
		return exeHome + "/" + binSubDir + "/" + exeName
	}
}
