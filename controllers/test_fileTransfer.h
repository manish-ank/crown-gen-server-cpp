#pragma once

#include <drogon/HttpController.h>
#include <drogon/HttpResponse.h>
using namespace drogon;

namespace test
{
  class fileTransfer : public drogon::HttpController<fileTransfer>
  {
  public:
    METHOD_LIST_BEGIN
    // use METHOD_ADD to add your custom processing function here;
    // METHOD_ADD(fileTransfer::get, "/{2}/{1}", Get); // path is /test/fileTransfer/{arg2}/{arg1}
    // METHOD_ADD(fileTransfer::your_method_name, "/{1}/{2}/list", Get); // path is /test/fileTransfer/{arg1}/{arg2}/list
    // ADD_METHOD_TO(fileTransfer::your_method_name, "/absolute/path/{1}/{2}/list", Get); // path is /absolute/path/{arg1}/{arg2}/list
    ADD_METHOD_TO(fileTransfer::file_transfer, "/curvature", Post, Options); // path is /absolute/path/{arg1}/{arg2}/list

    METHOD_LIST_END
    // your declaration of processing function maybe like this:
    // void get(const HttpRequestPtr& req, std::function<void (const HttpResponsePtr &)> &&callback, int p1, std::string p2);
    // void your_method_name(const HttpRequestPtr& req, std::function<void (const HttpResponsePtr &)> &&callback, double p1, int p2) const;
    void file_transfer(const HttpRequestPtr &req, std::function<void(const HttpResponsePtr &)> &&callback) const;
  };
}
