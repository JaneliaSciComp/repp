package repp

import (
	"os"
	"path/filepath"

	"go.uber.org/multierr"
	"golang.org/x/exp/maps"
)

func CollectFiles(locations []string) ([]string, error) {
	var allErrs error
	allFiles := map[string]string{}

	for _, location := range locations {
		files, err := collectFilesFromPathLocation(location)
		if err != nil {
			allErrs = multierr.Append(allErrs, err)
		} else {
			for _, f := range files {
				allFiles[f] = f
			}
		}
	}

	return maps.Values(allFiles), allErrs
}

func collectFilesFromPathLocation(path string) (files []string, err error) {
	if path != "" {
		pathInfo, statErr := os.Stat(path)
		if statErr != nil {
			// some error occurred - I don't care what type of error
			return files, statErr
		}
		if pathInfo.IsDir() {
			// collect files from sequence dir
			dirContent, dirErr := os.ReadDir(path)
			if dirErr != nil {
				return files, dirErr
			}
			for _, f := range dirContent {
				if !f.IsDir() {
					filePath, fpErr := filepath.Abs(filepath.Join(path, f.Name()))
					if fpErr != nil {
						return files, fpErr
					}
					files = append(files, filePath)
				}
			}
		} else if !pathInfo.IsDir() {
			filePath, fpErr := filepath.Abs(path)
			if fpErr != nil {
				return files, fpErr
			}
			files = append(files, filePath)
		}
	}
	return files, nil
}
