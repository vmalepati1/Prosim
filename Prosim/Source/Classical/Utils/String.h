#pragma once

#include <string>
#include <vector>

typedef std::string String;

namespace classical {

	namespace utils {

		template <typename ...Args>
		String StringWithFormat(const String& format, Args && ...args)
		{
			auto size = std::snprintf(nullptr, 0, format.c_str(), std::forward<Args>(args)...);
			String output(size, '\0');
			std::sprintf(&output[0], format.c_str(), std::forward<Args>(args)...);
			return output;
		}

		void CStringUpperCase(char *str);
		void CStringLowerCase(char *str);

		std::vector<String> SplitString(const String& string, const String& delimiters);
		std::vector<String> SplitString(const String& string, const char delimiter);
		std::vector<String> Tokenize(const String& string);
		std::vector<String> GetLines(const String& string);

		const char* FindToken(const char* str, const String& token);
		const char* FindToken(const String& string, const String& token);
		int FindStringPosition(const String& string, const String& search, unsigned int offset = 0);
		String StringRange(const String& string, unsigned int start, unsigned int length);
		String RemoveStringRange(const String& string, unsigned int start, unsigned int length);

		String GetBlock(const char* str, const char** outPosition = nullptr);
		String GetBlock(const String& string, unsigned int offset = 0);

		String GetStatement(const char* str, const char** outPosition = nullptr);

		bool StringContains(const String& string, const String& chars);
		bool StartsWith(const String& string, const String& start);
		int NextInt(const String& string);
		double ToDouble(const String &string);

	}

}