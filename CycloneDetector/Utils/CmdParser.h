#ifndef CMD_PARSER_H
#define CMD_PARSER_H

#include <string>
#include <variant>
#include <unordered_map>
#include <vector>
#include <type_traits>

//==========================================================================
//==========================================================================
//==========================================================================

class CmdParser 
{
public:
	using StringType = std::string;

protected:	
	using AvailableTypes = std::variant<
		bool *, 
		int *, double *, StringType *,
		std::vector<int> *, std::vector<double> *, std::vector<StringType> *
	>;

	/*
	template <typename T>
	using EnabledTypes = typename std::enable_if<
		std::is_same<T, int>::value ||
		std::is_same<T, double>::value ||
		std::is_same<T, std::string>::value
	>::type;
	*/

	template <typename T, typename... Ts>
	using EnabledTypes = typename std::enable_if<std::disjunction_v<std::is_same<T, Ts>... >>;
		
public:	
	struct ParamInfo 
	{
		ParamInfo() {}
		ParamInfo(const StringType & switchName, bool optional) : 
			switchName(switchName),
			desc(""),
			optional(optional)
		{}

		ParamInfo(const StringType & switchName, const StringType & desc, bool optional) :
			switchName(switchName),
			desc(desc),
			optional(optional)
		{}

		StringType switchName;
		StringType desc;
		bool optional;
	protected:
		int foundCount;		
		AvailableTypes value;

		friend class CmdParser;
	};

	static bool SaveArguments(const char * fileName, int argc, char ** argv);
	static CmdParser CreateFromSaved(const char * fileName);

	CmdParser(int argc, char ** argv);
	CmdParser(const CmdParser & p) = default;
	CmdParser(CmdParser && p);
	~CmdParser();

	CmdParser& operator=(const CmdParser&) = default;
	
	const std::vector<StringType> & GetOtherParams() const;
	void ClearTargets();

	template <typename T, typename = EnabledTypes<T, int, double, bool, std::string>>		
	void AddTarget(const char * switchName, T * target);
			
	template <typename T, typename = EnabledTypes<T, int, double, bool, std::string>>
	void AddTarget(const ParamInfo & info, T * target);
	

	template <typename T, typename = EnabledTypes<T, int, double, std::string>>
	void AddTarget(const char * switchName, std::vector<T> * target);

	template <typename T, typename = EnabledTypes<T, int, double, bool, std::string>>
	void AddTarget(const ParamInfo & info, std::vector<T> * target);

	void PrintHelp();
	
	bool Parse();

protected:
		
	int argc;
	char ** argv;
	bool freeArgv;
	bool allTargetsAreOptional;
	
	std::unordered_map<StringType, ParamInfo> params;
	std::vector<StringType> otherParams;

	CmdParser(int argc, char ** argv, bool needFreeArgv);

	template <typename T>
	T * GetVariadic(AvailableTypes & v);

};

//==========================================================================
//==========================================================================
//==========================================================================

template <typename T, typename>
void CmdParser::AddTarget(const char * switchName, T * target)
{
	ParamInfo pi;
	pi.desc = "";
	pi.optional = true;
	pi.switchName = switchName;
	this->AddTarget(pi, target);
};

template <typename T, typename>
void CmdParser::AddTarget(const ParamInfo & info, T * target)
{
	if constexpr (std::is_same<T, bool>::value)
	{
		//default is set to false
		*target = false;
	}
	auto ins = this->params.insert(std::pair(info.switchName, info)).first;
	ins->second.value = AvailableTypes(target);
	ins->second.foundCount = 0;

	this->allTargetsAreOptional &= info.optional;
};

template <typename T, typename>
void CmdParser::AddTarget(const char * switchName, std::vector<T> * target)
{
	ParamInfo pi;
	pi.desc = "";
	pi.optional = true;
	pi.switchName = switchName;
	this->AddTarget(pi, target);
};

template <typename T, typename>
void CmdParser::AddTarget(const ParamInfo & info, std::vector<T> * target)
{	
	auto ins = this->params.insert(std::pair(info.switchName, info)).first;
	ins->second.value = AvailableTypes(target);
	ins->second.foundCount = 0;

	this->allTargetsAreOptional &= info.optional;
};

/// <summary>
/// Helper method to get added target in a correct
/// data type from std::variant
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
template <typename T>
T * CmdParser::GetVariadic(AvailableTypes & v)
{
	if (std::holds_alternative<T *>(v))
	{
		return std::get<T *>(v);		
	}

	return nullptr;
}

#endif
