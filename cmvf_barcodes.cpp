#include <iostream>
#include <string>
#include "mvf_module/mvf_module.hpp"


int main(int argc, char** argv) {
  try {
    std::vector<std::string> paths;
    std::string order_path;
    bool positional = false;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (!positional && (arg == "--help" || arg == "-h")) { std::cout << "usage: cmvf_barcodes [--order simplex_order.txt] [--] field_0000.smp ...\n"; return 0; }
      if (!positional && arg == "--") { positional = true; continue; }
      if (!positional && arg == "--order") { if (!order_path.empty() || i + 1 == argc || !*argv[i + 1]) throw std::runtime_error("provide one simplex-order file after --order"); order_path = argv[++i]; continue; }
      if (!positional && !arg.empty() && arg.front() == '-') throw std::runtime_error("unknown option: " + arg);
      paths.push_back(arg);
    }
    if (paths.empty()) throw std::runtime_error("provide at least one .smp field; use --help for usage");
    cmvf::barcode_engine engine(cmvf::read_smp(paths.front()), order_path.empty() ? std::vector<int>{} : cmvf::read_order(order_path));
    for (std::size_t i = 1; i < paths.size(); ++i) {
      try { engine.advance(cmvf::read_smp(paths[i])); }
      catch (const std::exception& error) { throw std::runtime_error("field " + std::to_string(i) + " (" + paths[i] + "): " + error.what()); }
    }
    const auto bars = engine.barcodes();
    std::cout << "{\"schema\":\"mvf_barcodes_v1\",\"barcodes\":[";
    for (std::size_t i = 0; i < bars.size(); ++i) {
      const auto& bar = bars[i];
      if (i) std::cout << ',';
      std::cout << "{\"dimension\":" << bar.dimension << ",\"birth\":" << bar.birth << ",\"death\":";
      if (bar.death < 0) std::cout << "null"; else std::cout << bar.death;
      std::cout << '}';
    }
    std::cout << "]}\n";
    if (!std::cout) throw std::runtime_error("output write failed");
    return 0;
  } catch (const std::exception& error) { std::cerr << "cmvf_barcodes: " << error.what() << '\n'; return 1; }
}
