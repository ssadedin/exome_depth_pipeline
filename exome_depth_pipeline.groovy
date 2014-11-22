about title : "Simple exome depth pipeline that processes each chromosome in parallel"

load 'exome_depth.groovy'

run {
    chr(1..22,'X') * [ exome_depth ] + merge_ed
}
