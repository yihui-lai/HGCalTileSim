@class usr::ArgumentExtender
@ingroup Common
@details
Often in analysis programs, we need a short string for process indentification,
but a longer string for the publication/presentation outputs or change a
parameter in calculations calculations. An example of this case would be when
an analysis is split into muon/electron channels, where the calculations are
nearly entirely the same, but require the plotting to put in unique labels
for the two channel, while also potentially using slightly different numbers
in calculation (scale factors, integrated luminosity... etc).
Such behaviour could be listed in a json file as something like:
something like:

```json
{
  channel:{
    mu:{
      fullname:"#mu channel"
      plotname:"#mu+jets"
      lumi: 1050
    },
    el: {
      fullname:"#it{e} channel"
      plotname: "#it{e}+jets"
      lumi: 1000
    }
  }
}
```

Our class would need to first read the input arguments, and be able to get the
extended processes in a simple method. In the example above, we would be able
to call the program with the options "Options1":

```sh
./my_analysis_program --channel mu
```

And in the analysis code, the extended values could be obtained as:

```cpp
lumi     = args.ArgExt<double>("channel", "lumi"); // 1050
plotname = args.ArgExt<string>("channel","plotname"); // "@mu+jets"
arg     = args.Arg<string>("channel"); // "mu"
```

---

The ArgumentExtender class is a wrapper for:

- `boost::property_tree` to store the json file structure described above
  file.
- `boost::program_options::options_description` instance for defining the
  the the options to be parsed.
- `boost::program_options::variables_map` instance for storing the results
  of the program input options.

Functions are added to simplify the argument parsing process, as well as
getting the extended argument values. Additional parsing rules will be added to
help ensure that the json file formats between options are in sync.

In general, the typical usage of the class should be something like this:

1. Defining the class instance with json file(s). During this step, the parser
   also checks to see that the json file structure is consistent for each of the
   potential options (i.e. they have the same extended variables defined in the
   file.)
2. Defining the options_description for defining the options.
3. Take in the program inputs. During this step, the class will also check to
   see if the option is consistent with the json file (i.e. if the user input
   has its entry corresponding entry in the json file). The options defined
   in the json file will be assumed to be madatory argument.
4. Throughout the program, call the class to get the program inputs and
   extended variables.

For the use-case above. It would look something like:

```cpp
  usr::po::options_description desc("test options");
  desc.add_options()
    ("channel,c",usr::po::value<std::string>())
  ;

  usr::ArgumentExtender myargs("name1.json","name2.json");
  myargs.AddOptions(desc);
  myargs.ParseOptions(argc,argv);

  myargs.ArgExt<double>("channel","lumi");
```

---

## File path generation by options

In the typical case that a user running a command with various options:

```
my_analysis  --method Random --inputtype electron --fittype unbinned --minrange 0.02
```

Typically an analysis will want to generate files based on a myriad of input
type. So the target filename might look something like:

```
results/electron/MyPlot_MethodRandom_unbinned_min0p02.pdf
```

This class allows one to do this sort of file name specification via the
embedded `ArgPathScheme` class. To use this functionally one will typically need
to specify:

1. The prefix path of the generated file name
2. A list of options that will be used in the director of the file path
3. A list of options that will be used in the basename of the file path
4. The prefix of the basename.
5. The extension of the file

Using the above command as an example, the options description will be something
like:

```cpp
desc.add_options()
  ("method",value<std::string>() )
  ("inputtype", value<std::string>() )
  ("fittype", value<std::string>() )
  ("minrange", value<double>() )
;
```

The generation of the file name will need to be set up as something like:

```cpp
// Setting up a total directory range
arg.SetFilePrefix( "result/" )

arg.SetDirScheme( {"inputtype", ""} ); // Specifying input  to user for directory naming
                                       // The second string is left empty so just the
                                       // argument value will be used for path generation

arg.SetNameScheme( { // Specifying inputs to user for basename generation
  {"method", "Method"}, // Add prefix "Method" in front of argument value for path generation
  {"fittype", ""} ,     // No prefix to be added
  {"minrange", "min"}   // prefix "min" is to be added, double type will have slight differences.
});

auto f = arg.MakeFileName("MyPlot", "pdf"); // "MyPlot" is the prefix of the file name
                                            // with a pdf file extension.
auto f = arg.MakeFileName("MyOtherPlot", "pdf"); // Can now be used multiple times.
```

The scheme is effectively a list of two strings, the first being which option to
use, the second the string to be used in the final path name (can be blank).
This will generate the file name based on the input arguments and the defined
strings. Note that the raw, string-like input will be used, so if you input a
integer argument as as 1e6, and another time like 100000, this will generate
different file names, even if the setup is entirely the same. (Strictly
speaking, it will not actually be a list of two string, but for most use cases,
this will what it will look like.)

The directory scheme is to use different inputs to generated nested directory
names, while the name scheme will generate the resulting base name of the file.
A few things points to clean up is that decimal points will be replaced with the
character 'p', and list arguments will be separated by the '-' character. The
character separating the option inputs values will simply be '_'. If the
argument doesn't exists, the the string is simply not generated.

A bunch of functions is also available for generating file paths with commonly
used file extensions.
