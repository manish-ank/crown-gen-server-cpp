#include "test_fileTransfer.h"
#include <igl/readSTL.h>
#include "mesh.h"
#include <filesystem>
#include <fstream>
using namespace test;

std::string createZipFile(const std::vector<std::string> &filePaths, const std::string &outputZipPath)
{
    // In a real application, you'd use a library like minizip or another
    // to actually create the ZIP file.
    // For now, let's just create a dummy file to represent the zip.
    // This part requires external library integration.
    std::ofstream dummyZip(outputZipPath, std::ios::binary);
    if (dummyZip.is_open())
    {
        dummyZip << "This is a dummy zip file containing:\n";
        for (const auto &path : filePaths)
        {
            dummyZip << "- " << path << "\n";
            // In a real zip, you'd add the file content here
            // (e.g., read from 'path' and write to 'dummyZip')
        }
        dummyZip.close();
        return outputZipPath;
    }
    return ""; // Return empty string on failure
}
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
    std::string file_extension = fileParser.getFiles()[0].getFileExtension().data();
    fileParser.getFiles()[0].save("./");
    fileParser.getFiles()[1].save("./");

    std::string file_1_name = fileParser.getFiles()[0].getFileName();
    std::string file_2_name = fileParser.getFiles()[1].getFileName();
    // ############################################################3###########################################################################################
    CrownGen::Mesh mesh1, mesh2;
    std::string out_file_1_name = "data1.obj";
    std::string out_file_2_name = "data2.obj";
    std::string color_1_json = "data_1_colors.json";
    std::string color_2_json = "data_2_colors.json";

    LOG_INFO << "init";
    mesh1.init(file_1_name, out_file_1_name,color_1_json);
    mesh2.init(file_2_name, out_file_2_name,color_2_json);
    LOG_INFO << "curvcalc";

    mesh1.curvatureCalculation(0.0);
    mesh2.curvatureCalculation(0.0);
    LOG_INFO << "adjcalc";

    mesh1.calculateAdjacencyVerticesFaces();
    mesh2.calculateAdjacencyVerticesFaces();

    LOG_INFO << "seg";
    mesh1.regionGrowingSegmentation(30, 0.5, 10, 100, 100000);
    mesh2.regionGrowingSegmentation(30, 0.5, 10, 100, 100000);

    LOG_INFO << "zip";

    // ############################################################3###########################################################################################
    // zip the file and send it as response
    std::string zipname = "./files.zip";

    if (!std::filesystem::exists(out_file_1_name) || !std::filesystem::is_regular_file(out_file_1_name) || !std::filesystem::exists(out_file_2_name) || !std::filesystem::is_regular_file(out_file_2_name) ||
        !std::filesystem::exists(color_1_json) || !std::filesystem::is_regular_file(color_1_json) ||
        !std::filesystem::exists(color_2_json) || !std::filesystem::is_regular_file(color_2_json))
    {
        auto resp = HttpResponse::newHttpResponse();
        resp->setStatusCode(k404NotFound);
        resp->setBody("One or multiple source files not found.");
        callback(resp);
        return;
    }

    std::vector<std::string> filesToZip = {out_file_1_name, out_file_2_name, color_1_json, color_2_json};
    std::string generatedZipPath = createZipFile(filesToZip, zipname);

    if (generatedZipPath.empty())
    {
        auto resp = HttpResponse::newHttpResponse();
        resp->setStatusCode(k500InternalServerError);
        resp->setBody("Failed to create ZIP archive.");
        callback(resp);
        return;
    }

    // 3. Send the ZIP file as a response
    auto resp = HttpResponse::newFileResponse(generatedZipPath,
                                              0,
                                              0,
                                              true,
                                              "segmented.zip",             // Custom name for the downloaded zip
                                              drogon::CT_APPLICATION_ZIP); // Set content type to zip

    // ############################################################3###########################################################################################
    callback(resp);
}
