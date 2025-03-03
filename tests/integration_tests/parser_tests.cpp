#include "gtest/gtest.h"

#include "src/input/Parser.h"

#include <filesystem>

std::vector<std::filesystem::path> list_of_files = {
	"10.1021_ja0305959/d-nit2.yml",
	"10.1039_B414284E/3d-nit2.yml",
	"10.1039_B414284E/dtd-nit2.yml",
	"10.1039_D0DT03184D/manganese_cocrystal.yml",
	"10.1039_D0DT03184D/nickel_crystal_first_model.yml",
	"10.1039_D0DT03184D/nickel_crystal_second_model.yml",
	"10.1039_D0DT03184D/nickel_crystal_third_model.yml"
};

TEST(integration_parser_tests, do_not_throw_on_examples) {
	const auto old_path = std::filesystem::current_path();
	common::Logger::set_level(common::PrintLevel::off);
	for (const std::string& file : list_of_files) {
		// this "../../" may cause troubles in the case of unusual location of build directory 
		auto absolute_path = old_path.parent_path().parent_path() / "examples" / file;
		std::filesystem::current_path(absolute_path.parent_path());
		EXPECT_NO_THROW(auto parser = input::Parser(absolute_path));
		common::Logger::set_level(common::PrintLevel::off);
	}
	std::filesystem::current_path(old_path);
}