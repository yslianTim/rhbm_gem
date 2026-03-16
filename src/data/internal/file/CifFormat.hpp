#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class ModelObject;
struct CifFormatState;

class CifFormat {
  public:
    CifFormat();
    ~CifFormat();
    std::unique_ptr<ModelObject> ReadModel(std::istream& stream, const std::string& source_name);
    void WriteModel(const ModelObject& model_object, std::ostream& stream, int par);

  private:
    std::unique_ptr<CifFormatState> m_read_state;
};

} // namespace rhbm_gem
