
import subprocess
# import random


flagstat = "00.flagstat.sh"
filter_cram = "01.filter_cram.sh"
subset_region = "02.subset_cram_by_region.sh"
split_by_read_orientation = "03.split_by_read_orientation.sh"
merge_read_orientations = "03b.merge_read_orientations.sh"
pileup = "04.pileup.sh"


def submit_workflow(input_cram, reference, bed_forward, bed_reverse):
    # Step 0: Flagstat on input CRAM.
    submit_job(flagstat, i=input_cram)

    # Step 1: Filter input CRAM.
    step_1 = submit_job(filter_cram,
                        r=reference,
                        i=input_cram,
                        o=(filtered := input_cram.replace(".cram", ".filtered.cram")))
    submit_job(flagstat, i=filtered, dependency_ids=step_1)

    # Step 2: Subset CRAM by transcription regions.
    step_2 = dict()
    for region, region_name in zip([bed_forward, bed_reverse],
                                   ["forward_coding", "reverse_coding"]):
        out = filtered.replace(".cram", f".{region_name}.cram")
        substep = submit_job(subset_region,
                             dependency_ids=step_1,
                             r=reference,
                             i=filtered,
                             b=region,
                             o=out)
        submit_job(flagstat, i=out, dependency_ids=substep)
        step_2[substep] = out

    # Step 3: Split CRAMS by read orientation.
    step_3 = []
    orientations = ["F1", "R2", "F2", "R1"]

    for substep2_id, cram in step_2.items():
        orients = dict()
        for orientation in orientations:
            out = cram.replace(".cram", f".{orientation}.cram")
            substep = submit_job(split_by_read_orientation,
                                 dependency_ids=substep2_id,
                                 r=reference,
                                 i=cram,
                                 o=out,
                                 d=orientation)
            orients[orientation] = (substep, out)
            submit_job(flagstat, i=out, dependency_ids=substep)
        step_3.append(orients)

    # Step 3b: Step Merge F1R2s and F2R1s.
    step_3b = dict()
    for strand_dict in step_3:
        for (a, b) in [orientations[:2], orientations[2:]]:
            a_step, a_cram = strand_dict[a]
            b_step, b_cram = strand_dict[b]
            out = ".".join(a_cram.split(".")[:-2]) + f".{a}{b}.cram"
            substep = submit_job(merge_read_orientations,
                                 dependency_ids=[a_step, b_step],
                                 a=a_cram,
                                 b=b_cram,
                                 o=out)
            step_3b[substep] = out
            submit_job(flagstat, i=out, dependency_ids=substep)

    step_3 = {substep: out for orients in step_3 for _, (substep, out) in orients.items()}
    step_3.update(step_3b)

    # Step 4: Produce pileups.
    for substep, cram in step_3.items():
        submit_job(pileup,
                   dependency_ids=substep,
                   r=reference,
                   i=cram)


def submit_job(tool, dependency_ids=None, dependency_mode="afterok", **kwargs):
    dependencies = parse_dependencies(dependency_ids, dependency_mode)
    submission = ["sbatch", "--parsable", dependencies] + make_command(tool, **kwargs)
    submission = [piece for piece in submission if piece]
    submission = subprocess.run(submission, capture_output=True, text=True, check=True)
    jobid = int(submission.stdout.strip())
    # jobid = random.randint(0, 1000)
    # print(f"{tool}\t{jobid}\t{' '.join(submission)}")
    return jobid


def parse_dependencies(dependency_ids=None, dependencies_mode="afterok"):
    if dependency_ids is None:
        return None
    if isinstance(dependency_ids, int):
        dependency_ids = [dependency_ids]
    dependencies = [str(dependency_id) for dependency_id in dependency_ids]
    dependencies = f"--dependency={dependencies_mode}:" + ",".join(dependencies)
    return dependencies


def make_command(tool, **kwargs):
    flags = [f"-{key}" for key in kwargs.keys()]
    options = [f"{arg}" for arg in kwargs.values()]
    command = [tool] + [item for pair in zip(flags, options) for item in pair]
    return command


if __name__ == '__main__':
    import sys
    submit_workflow(*sys.argv[1:])
