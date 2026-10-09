#!/bin/bash
source /home/chizixin/OpenFOAM/OpenFOAM-10/etc/bashrc
# 汇总日志
SUMMARY="summary_OptRodas34_A8R1_Mass_Fast.log"
echo "Summary generated at $(date)" > "$SUMMARY"
echo "-------------------------------------" >> "$SUMMARY"
# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions
# 遍历所有包含 system/controlDict 的算例目录
for case in */ ; do
    if [[ -f "${case}/system/controlDict" ]]; then
        echo "Processing case: $case"

        cd "$case"
        foamDictionary -entry "jacobian" -set "fast" constant/chemistryProperties >/dev/null
        foamDictionary -entry "massFraction" -set "true" constant/chemistryProperties >/dev/null
        foamDictionary -entry "chemistryType/solver" -set "OptRodas34" constant/chemistryProperties >/dev/null
        foamDictionary -entry "chemistryType/method" -set "FastChemistryModel" constant/chemistryProperties >/dev/null
        foamDictionary -entry "OptRodas34Coeffs/absTol" -set "1e-8" constant/chemistryProperties >/dev/null
        foamDictionary -entry "OptRodas34Coeffs/relTol" -set "1e-1" constant/chemistryProperties >/dev/null
        # 运行 chemFoam 并生成 log
        chemFoam > log 2>&1

        # 提取倒数第三行
        last3=$(tail -n 3 log | head -n 1)

        cd ..

        # 写入 summary.log
        echo "$case : $last3" >> "$SUMMARY"
    fi
done

echo "Done. Summary written to $SUMMARY"

