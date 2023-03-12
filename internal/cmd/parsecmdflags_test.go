package cmd

import (
	"reflect"
	"testing"

	"github.com/spf13/cobra"
)

func createTestCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "sequence",
		Short: "Test sequence command",
	}
	cmd.Flags().String("exclude", "", "filters")
	cmd.Flags().String("dbs", "", "dbnames")
	return cmd
}

func Test_getFilters(t *testing.T) {
	cmd := createTestCmd()
	tests := []struct {
		name       string
		filtersArg string
		want       []string
	}{
		{
			"biobrick separated from year by commas",
			"tests,BBa_k222000   ,  2004",
			[]string{"TESTS", "BBA_K222000", "2004"},
		},
		{
			"single year",
			"2004",
			[]string{"2004"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			cmd.SetArgs([]string{
				"sequence",
				"--exclude",
				tt.filtersArg,
			})
			cmd.Run = func(cmd *cobra.Command, args []string) {
				got := extractExcludedValues(cmd)
				if !reflect.DeepEqual(got, tt.want) {
					t.Errorf("parse filters = %v, want %v", got, tt.want)
				}
			}
			err := cmd.Execute()
			if err != nil {
				t.Fail()
			}
		})
	}
}

func Test_getDBs(t *testing.T) {
	cmd := createTestCmd()
	tests := []struct {
		name       string
		filtersArg string
		want       []string
	}{
		{
			"comma separated db names",
			"db1,  Db2   ,  dB3,DB4",
			[]string{"db1", "Db2", "dB3", "DB4"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			cmd.SetArgs([]string{
				"sequence",
				"--dbs",
				tt.filtersArg,
			})
			cmd.Run = func(cmd *cobra.Command, args []string) {
				got := extractDbNames(cmd)
				if !reflect.DeepEqual(got, tt.want) {
					t.Errorf("parse dbs = %v, want %v", got, tt.want)
				}
			}
			err := cmd.Execute()
			if err != nil {
				t.Fail()
			}
		})
	}
}

func Test_guessOutput(t *testing.T) {
	type args struct {
		in string
	}
	tests := []struct {
		name    string
		args    args
		wantOut string
	}{
		{
			"append json suffix",
			args{
				in: "./test_file.fa",
			},
			"./test_file.output.json",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotOut := guessOutput(tt.args.in); gotOut != tt.wantOut {
				t.Errorf("guessOutput() = %v, want %v", gotOut, tt.wantOut)
			}
		})
	}
}

func Test_combineFeaturesIntoCSV(t *testing.T) {
	tests := []struct {
		name    string
		args    []string
		wantOut string
	}{
		{
			"combine comma separated features",
			[]string{
				"p10 promoter, mEGFP",
				"T7 terminator",
			},
			"p10 promoter,mEGFP T7 terminator",
		},
		{
			"combine space separated features",
			[]string{
				"p10   promoter",
				"mEGFP T7    terminator",
			},
			"p10   promoter mEGFP T7    terminator",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotOut := combineAllIntoCSV(tt.args); gotOut != tt.wantOut {
				t.Errorf("combineAllIntoCSV() = %v, want %v", gotOut, tt.wantOut)
			}
		})
	}

}
