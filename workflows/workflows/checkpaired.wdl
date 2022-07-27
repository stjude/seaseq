version 1.0
import "../tasks/sratoolkit.wdl" as sra

workflow checkpaired {

    input {
        Array[String] sraid
    }

    scatter (eachsra in sraid) {
        call sra.srameta {
            input :
               sra_id=eachsra
        }
    }

    output {
        Array[Boolean] paired = srameta.paired_end
    }
}
