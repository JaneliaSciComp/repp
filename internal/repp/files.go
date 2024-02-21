package repp

import (
	"os"
	"path/filepath"
	"slices"
	"strings"

	"go.uber.org/multierr"
	"golang.org/x/exp/maps"
)

func CollectFiles(locations []string) ([]string, error) {
	var allErrs error
	var expandedLocations []string

	for _, l := range locations {
		if strings.TrimSpace(l) == "" {
			continue
		}
		paths, err := filepath.Glob(l)
		if err != nil {
			rlog.Infof("Error expanding %s: %v", l, err)
			continue
		}
		expandedLocations = append(expandedLocations, paths...)
	}

	allFiles := map[string]string{}
	for _, l := range expandedLocations {
		files, err := collectFilesFromPathLocation(l)
		if err != nil {
			allErrs = multierr.Append(allErrs, err)
		} else {
			for _, f := range files {
				allFiles[f] = f
			}
		}
	}

	// sort the final results
	allFilePaths := maps.Values(allFiles)
	fnameCmp := func(fp1, fp2 string) int {
		fn1 := filepath.Base(fp1)
		fn2 := filepath.Base(fp2)
		return strings.Compare(fn1, fn2)
	}
	slices.SortFunc(allFilePaths, fnameCmp)

	return allFilePaths, allErrs
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
