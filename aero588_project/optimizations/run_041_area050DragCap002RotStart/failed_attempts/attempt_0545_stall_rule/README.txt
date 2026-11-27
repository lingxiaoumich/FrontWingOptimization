run_041's first attempt (05:45) was ended after 20 evaluations by run_chain.py's stall rule: its start is 76 % over
    the drag cap, PSQP traded downforce for drag, and the rule saw no C_L c gain in 15 min. run_chain.py now counts only
    iterates within 1 % of the caps (user, 2026-11-25: "fix rule, rerun"). The attempt is kept in failed_attempts/; the
    rerun takes lane A once this queue's first run (run_040) has left it.
