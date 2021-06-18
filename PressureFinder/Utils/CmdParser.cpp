#include "./CmdParser.h"

#include "../Macros.h"

//==========================================================================
//==========================================================================
//==========================================================================

/// <summary>
/// Save command line arguments to file
/// This file can be loaded to parser
/// Each argument - separate line
/// </summary>
/// <param name="fileName"></param>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
bool CmdParser::SaveArguments(const char * fileName, int argc, char ** argv)
{
	FILE * fp;
	my_fopen(&fp, fileName, "w");
	if (fp == nullptr)
	{		
		return false;
	}

	int i = 1;	
	while (i < argc)
	{
		fprintf(fp, "%s\n", argv[i]);		
		i++;
	}

	fclose(fp);
	return true;
}

/// <summary>
/// Load command line arguments saved in file
/// Each argument - separate line 
/// Max line length 1024 bytes
/// </summary>
/// <param name="fileName"></param>
/// <returns></returns>
CmdParser CmdParser::CreateFromSaved(const char * fileName)
{
	FILE * fp;
	my_fopen(&fp, fileName, "r");
	if (fp == nullptr)
	{
		return CmdParser(0, nullptr);
	}

	int argc = 1;
	char c;
	for (c = getc(fp); c != EOF; c = getc(fp))
	{
		if (c == '\n')
		{
			argc++;
		}
	}

	fseek(fp, 0, SEEK_SET);

	char ** argv = new char *[argc + 1];
	argv[argc] = nullptr;

	char * buff = new char[1024];
	strcpy(buff, "from file");
	argv[0] = buff;

	for (int i = 1; i <= argc; i++)
	{
		buff = new char[1024];
		fgets(buff, 1024, (FILE*)fp);

		int len = static_cast<int>(strlen(buff));
		if (len > 0 && buff[len - 1] == '\n')
		{
			buff[len - 1] = 0;
			len--;
		}
		if (len > 0 && buff[len - 1] == '\r')
		{
			buff[len - 1] = 0;
		}

		argv[i] = buff;
	}

	fclose(fp);

	return CmdParser(argc, argv, true);
}

//==========================================================================
//==========================================================================
//==========================================================================


/// <summary>
/// ctor
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
CmdParser::CmdParser(int argc, char ** argv) : 
	CmdParser(argc, argv, false)
{
}

/// <summary>
/// ctor
/// params argv are "user" allocated and need to be released
/// after CmdParser finishes
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <param name="needFreeArgv"></param>
CmdParser::CmdParser(int argc, char ** argv, bool needFreeArgv) : 
	argc(argc),
	argv(argv),
	freeArgv(needFreeArgv),
	allTargetsAreOptional(true)
{
}

/// <summary>
/// move ctor
/// </summary>
/// <param name="p"></param>
CmdParser::CmdParser(CmdParser && p) :
	argc(p.argc),
	argv(p.argv),
	freeArgv(p.freeArgv),
	allTargetsAreOptional(p.allTargetsAreOptional)
{
	p.argv = nullptr;
	p.argc = -1;
}

/// <summary>
/// dtor
/// </summary>
CmdParser::~CmdParser()
{
	if ((freeArgv) && (argv != nullptr))
	{
		for (int i = 0; i <= argc; i++)
		{
			SAFE_DELETE_ARRAY(argv[i]);
		}

		SAFE_DELETE_ARRAY(argv);
	}
}

/// <summary>
/// Get parameters that were not assigned to be parsed directly
/// </summary>
/// <returns></returns>
const std::vector<CmdParser::StringType> & CmdParser::GetOtherParams() const
{
	return this->otherParams;
}

/// <summary>
/// Clear set up parameters
/// - you have to call AddTarget again before calling Parse
/// - if not, all params will be loaded to otherParams
/// </summary>
void CmdParser::ClearTargets()
{
	this->params.clear();
	this->otherParams.clear();
}

/// <summary>
/// Print help based on settings
/// </summary>
void CmdParser::PrintHelp()
{		

	printf("=================\n");
	printf("===== HELP ======\n");
	printf("=================\n");

	for (auto & p : this->params)
	{
		printf("-%s", p.second.switchName.c_str());

		if (std::holds_alternative<bool *>(p.second.value))
		{
			printf(" (default: false) ");
		}
		else if ((std::holds_alternative<int *>(p.second.value)) ||
			(std::holds_alternative<std::vector<int> *>(p.second.value)))
		{
			printf(" <int> ");
			if (std::holds_alternative<std::vector<int> *>(p.second.value))
			{
				printf(" (multivalue) ");
			}
		}
		else if ((std::holds_alternative<double *>(p.second.value)) ||
			(std::holds_alternative<std::vector<double> *>(p.second.value)))
		{
			printf(" <double> ");
			if (std::holds_alternative<std::vector<double> *>(p.second.value))
			{
				printf(" (multivalue) ");
			}
		}
		else if ((std::holds_alternative<StringType *>(p.second.value)) ||
			(std::holds_alternative<std::vector<StringType> *>(p.second.value)))
		{
			printf(" <string> ");

			if (std::holds_alternative<std::vector<StringType> *>(p.second.value))
			{
				printf(" (multivalue) ");
			}
		}

		
		if (p.second.foundCount != 0)
		{			
			printf(" (current value: ");

			if (GetVariadic<bool>(p.second.value) != nullptr)
			{
				printf("true");
			}			
			else if (auto tmp = GetVariadic<int>(p.second.value))
			{
				printf("%d", *tmp);
			}
			else if (auto tmp = GetVariadic<double>(p.second.value))
			{
				printf("%f", *tmp);
			}
			else if (auto tmp = GetVariadic<StringType>(p.second.value))
			{
				printf("%s", tmp->c_str());
			}
			else 
			{
				for (int i = 0; i < p.second.foundCount; i++)
				{
					if (auto tmp = GetVariadic<std::vector<int>>(p.second.value))
					{
						printf("%d", tmp->at(i));
					}
					else if (auto tmp = GetVariadic<std::vector<double>>(p.second.value))
					{
						printf("%f", tmp->at(i));
					}
					else if (auto tmp = GetVariadic<std::vector<StringType>>(p.second.value))
					{
						printf("%s", tmp->at(i).c_str());
					}
					
					if (i + 1 != p.second.foundCount)
					{
						printf(", ");
					}
				}
			}
			

			printf(") ");
		}

		if (p.second.optional)
		{
			printf(" (optional) ");
		}
		else if (p.second.foundCount == 0)
		{
			printf(" (!!! parametr missing !!!) ");
		}

		if (p.second.desc != "")
		{
			printf("\n\tDescription: %s", p.second.desc.c_str());
		}
		printf("\n\n");
	}
	
	if (this->otherParams.size() != 0)
	{
		printf("=================\n");
		printf("Other unrecognized parameters:\n");
		for (auto & v : this->otherParams)
		{
			printf("%s\n", v.c_str());
		}
	}

	printf("=================\n");
}

/// <summary>
/// Parse command line arguments
/// If number of required parametrs is not met, returns false
/// </summary>
/// <returns></returns>
bool CmdParser::Parse()
{	
	if (argc == 0)
	{		
		return this->allTargetsAreOptional;
	}

	int i = 1;
	
	for (auto & p : this->params) 
	{
		p.second.foundCount = 0;
	}

	while (i < argc)
	{
		//first test without "-"
		auto it = this->params.find(argv[i]);
		if (it == this->params.end())
		{
			//not found - test with "-"
			it = this->params.find(argv[i] + 1);
			if (it == this->params.end())
			{
				this->otherParams.push_back(argv[i]);
				i++;
				continue;
			}			
		}
		
		it->second.foundCount++;
		i++;
		
		if (auto tmp = GetVariadic<bool>(it->second.value))
		{
			*tmp = true;
		}		
		else
		{
			const char * val = argv[i++];
			if (auto tmp = GetVariadic<int>(it->second.value))
			{
				*tmp = atoi(val);
			}
			else if (auto tmp = GetVariadic<double>(it->second.value))
			{
				*tmp = atof(val);
			}
			else if (auto tmp = GetVariadic<StringType>(it->second.value))
			{
				*tmp = val;
			}
			else if (auto tmp = GetVariadic<std::vector<int>>(it->second.value))
			{
				tmp->push_back(atoi(val));
			}
			else if (auto tmp = GetVariadic<std::vector<double>>(it->second.value))
			{
				tmp->push_back(atof(val));
			}
			else if (auto tmp = GetVariadic<std::vector<StringType>>(it->second.value))
			{
				tmp->push_back(val);
			}
		}
	}


	for (auto & p : this->params)
	{
		if ((p.second.optional == false) && (p.second.foundCount == 0))
		{
			return false;
		}
	}

	return true;
}