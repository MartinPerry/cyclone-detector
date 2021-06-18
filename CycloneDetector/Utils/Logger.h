#ifndef LOGGER_H
#define LOGGER_H


#ifndef MY_LOG_ERROR
#   define MY_LOG_ERROR(...) MyUtils::Logger::LogError(__VA_ARGS__)
#endif

#ifndef MY_LOG_WARNING
#   define MY_LOG_WARNING(...) MyUtils::Logger::LogWarning(__VA_ARGS__)
#endif

#ifndef MY_LOG_INFO
#	if defined(_DEBUG) || defined(DEBUG) || defined (USE_RELEASE_LOG_INFO)
#		define MY_LOG_INFO(...) MyUtils::Logger::LogInfo(__VA_ARGS__)
#	else
#		define MY_LOG_INFO(...) 
#	endif
#endif

#ifndef MY_LOG
#   define MY_LOG(...) MyUtils::Logger::LogMessage(__VA_ARGS__)
#endif

#include <cstdio>

namespace MyUtils
{
	class Logger
	{
		public:
			typedef enum LOG_OUTPUT { LOG_STDOUT = 0, LOG_FILE = 1, LOG_STDERR = 2 } LOG_OUTPUT;

			static void Initialize();
			static void Destroy();
			static Logger * GetInstance();

			static void LogError(const char * message, ...);
			static void LogWarning(const char * message, ...);
			static void LogInfo(const char * message, ...);
			static void LogMessage(const char * message, ...);


			void LogToFile(const char * fileName);
			void LogToStdout();
			void LogToStder();						

			void EnableErrorLogging(LOG_OUTPUT type);
			void EnableWarningLogging(LOG_OUTPUT type);
			void EnableInfoLogging(LOG_OUTPUT type);

			void DisableErrorLogging(LOG_OUTPUT type);
			void DisableWarningLogging(LOG_OUTPUT type);
			void DisableInfoLogging(LOG_OUTPUT type);

			bool IsErrorLoggingEnabled(LOG_OUTPUT type);
			bool IsWarningLoggingEnabled(LOG_OUTPUT type);
			bool IsInfoLoggingEnabled(LOG_OUTPUT type);

		private:
#ifdef _WIN32
			static const int COLOR_RED = 12;
			static const int COLOR_GREEN = 10;
			static const int COLOR_YELLOW = 14;
			static const int COLOR_WHITE = 15;
			static const int COLOR_BLUE = 9;
#else
            static const int COLOR_RED = 12;
            static const int COLOR_GREEN = 10;
            static const int COLOR_YELLOW = 14;
            static const int COLOR_WHITE = 15;
			static const int COLOR_BLUE = 9;
#endif
        
        
			static Logger * instanceLogger;

			Logger();
			~Logger();

			bool enableErrors[3];
			bool enableWarnings[3];
			bool enableInfo[3];

			FILE * loggerOutput[3];

			void StartColor(int colorID);
			void EndColor(int colorID);
			
			
	};
}


#endif
