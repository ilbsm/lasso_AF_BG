reinitialize
import sys

python
template = sys.argv[1]
model = sys.argv[2]

template_name = template.split("/")[-1].split(".")[0]
model_name = model.split("/")[-1].split(".")[0]
python end

cmd.load("%s"%template)
cmd.load("%s"%model)

alignment = cmd.super("%s"%template_name, "%s"%model_name)

python
print(str(alignment[0]))
with open("result_file.tsv", "a") as newfile:
    newfile.write("\t".join([template_name, model_name, str(alignment[0])]))
    newfile.write("\n")
python end


