#include "String.h"

namespace classical {

	namespace utils {

		void CStringUpperCase(char *str) {
			int c = 0;

			while (str[c] != '\0') {
				if (str[c] >= 'a' && str[c] <= 'z') {
					str[c] = str[c] - 32;
				}
				c++;
			}
		}

		void CStringLowerCase(char *str) {
			int c = 0;

			while (str[c] != '\0') {
				if (str[c] >= 'A' && str[c] <= 'Z') {
					str[c] = str[c] + 32;
				}
				c++;
			}
		}

		std::vector<String> SplitString(const String& string, const String& delimiters)
		{
			size_t start = 0;
			size_t end = string.find_first_of(delimiters);

			std::vector<String> result;

			while (end <= String::npos)
			{
				String token = string.substr(start, end - start);
				if (!token.empty())
					result.push_back(token);

				if (end == String::npos)
					break;

				start = end + 1;
				end = string.find_first_of(delimiters, start);
			}

			return result;
		}

		std::vector<String> SplitString(const String& string, const char delimiter)
		{
			return SplitString(string, String(1, delimiter));
		}

		std::vector<String> Tokenize(const String& string)
		{
			return SplitString(string, " \t\n");
		}

		std::vector<String> GetLines(const String& string)
		{
			return SplitString(string, "\n");
		}

		const char* FindToken(const char* str, const String& token)
		{
			const char* t = str;
			while (t = strstr(t, token.c_str()))
			{
				bool left = str == t || isspace(t[-1]);
				bool right = !t[token.size()] || isspace(t[token.size()]);
				if (left && right)
					return t;

				t += token.size();
			}
			return nullptr;
		}

		const char* FindToken(const String& string, const String& token)
		{
			return FindToken(string.c_str(), token);
		}

		int FindStringPosition(const String& string, const String& search, unsigned int offset)
		{
			const char* str = string.c_str() + offset;
			const char* found = strstr(str, search.c_str());
			if (found == nullptr)
				return -1;
			return (int)(found - str) + offset;
		}

		String StringRange(const String& string, unsigned int start, unsigned int length)
		{
			return string.substr(start, length);
		}

		String RemoveStringRange(const String& string, unsigned int start, unsigned int length)
		{
			String result = string;
			return result.erase(start, length);
		}

		String GetBlock(const char* str, const char** outPosition)
		{
			const char* end = strstr(str, "}");
			if (!end)
				return String(str);

			if (outPosition)
				*outPosition = end;
			unsigned int length = end - str + 1;
			return String(str, length);
		}

		String GetBlock(const String& string, unsigned int offset)
		{
			const char* str = string.c_str() + offset;
			return GetBlock(str);
		}

		String GetStatement(const char* str, const char** outPosition)
		{
			const char* end = strstr(str, ";");
			if (!end)
				return String(str);

			if (outPosition)
				*outPosition = end;
			unsigned int length = end - str + 1;
			return String(str, length);
		}

		bool StringContains(const String& string, const String& chars)
		{
			return string.find(chars) != String::npos;
		}

		bool StartsWith(const String& string, const String& start)
		{
			return string.find(start) == 0;
		}

		int NextInt(const String& string)
		{
			const char* str = string.c_str();
			for (unsigned int i = 0; i < string.size(); i++)
			{
				if (isdigit(string[i]))
					return atoi(&string[i]);
			}
			return -1;
		}

		double ToDouble(const String &string) {
			return atof(string.c_str());
		}

	}

}