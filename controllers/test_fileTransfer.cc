#include "test_fileTransfer.h"
#include <igl/readSTL.h>
#include "mesh.h"
using namespace test;

// Add definition of your processing function here

void fileTransfer::file_transfer(const HttpRequestPtr &req, std::function<void(const HttpResponsePtr &)> &&callback) const
{

    if (req->method() == drogon::HttpMethod::Options)
    {
        auto resp = HttpResponse::newHttpResponse();
        resp->setStatusCode(k200OK);
        resp->addHeader("Access-Control-Allow-Origin", "*");
        resp->addHeader("Access-Control-Allow-Methods", "POST, OPTIONS");
        resp->addHeader("Access-Control-Allow-Headers", "Content-Type");
        callback(resp);
        return;
    }

    MultiPartParser fileParser;
    fileParser.parse(req);

    if (fileParser.getFiles().empty())
    {
        LOG_INFO << "File Not found";
        auto resp = HttpResponse::newHttpJsonResponse(Json::Value("File not found\n"));
        resp->setStatusCode(k400BadRequest);
        callback(resp);
        return;
    }
    size_t num_of_files = fileParser.getFiles().size();
    LOG_INFO << "Number of files: " << num_of_files;

    fileParser.getFiles()[0].saveAs("./data.stl");

    //############################################################3###########################################################################################

    CrownGen::Mesh mesh;
    mesh.initSTL("data.stl");
    
    //############################################################3###########################################################################################


    auto resp = HttpResponse::newHttpResponse();
    resp->setStatusCode(k200OK);
    resp->setBody("File Recieved\n");
    resp->addHeader("Access-Control-Allow-Origin", "*");
    resp->addHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
    resp->addHeader("Access-Control-Allow-Headers", "Content-Type");
    callback(resp);
}
