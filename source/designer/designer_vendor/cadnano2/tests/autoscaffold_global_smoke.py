"""Manual regression runner for the global AutoCS scaffold route."""

import json
import faulthandler
import os
import sys

import cadnano2.cadnano as cadnano


def main(path):
    faulthandler.enable()
    faulthandler.dump_traceback_later(120, repeat=True)
    from PyQt6.QtWidgets import QApplication
    qtApp = QApplication.instance() or QApplication([])
    app = cadnano.initAppWithoutGui([])
    app.prefs.honeycombRows = 50
    app.prefs.honeycombCols = 50
    app.prefs.honeycombSteps = 2
    app.prefs.squareRows = 50
    app.prefs.squareCols = 50
    app.prefs.squareSteps = 2
    class DummySignal:
        def emit(self, *unused_args):
            pass
    app.documentWasCreatedSignal = DummySignal()
    from cadnano2.model.document import Document
    from cadnano2.model.io.decoder import decode
    document = Document()
    if path in ('--synthetic-triangle', '--synthetic-honeycomb-flat',
                '--synthetic-square-flat', '--synthetic-square-short54',
                '--synthetic-square-short32',
                '--synthetic-square-four',
                '--synthetic-square-grid128',
                '--synthetic-square-modules32'):
        lattice = ('square' if path.startswith('--synthetic-square')
                   else 'honeycomb')
        if hasattr(document, 'addGuidedPart'):
            part = document.addGuidedPart(
                lattice, 50 if path == '--synthetic-square-grid128' else 6,
                50 if path == '--synthetic-square-grid128' else 14,
                4 if path == '--synthetic-square-grid128' else
                4 if path == '--synthetic-square-four' else
                (9 if lattice == 'square' else 12))
        else:
            app.prefs.squareSteps = 4
            app.prefs.squareRows = 50
            app.prefs.squareCols = 50
            part = (document.addSquarePart() if lattice == 'square'
                    else document.addHoneycombPart())
        shortLength = (54 if path.endswith('short54') else
                       32 if path.endswith('short32') else None)
        records = ([
                    (40 + row, 30 + column, 0, 127)
                    for row in range(4) for column in range(8)]
                   if path == '--synthetic-square-modules32' else
                   [(47, 27, 0, 127), (47, 28, 0, 127),
                    (47, 29, 0, 127), (47, 30, 0, 127),
                    (48, 30, 0, 127), (48, 29, 0, 127),
                    (48, 28, 0, 127), (48, 27, 0, 127)]
                   if path == '--synthetic-square-grid128' else
                   [(1, column, 9, 119)
                    for column in range(2, 6)]
                   if path == '--synthetic-square-four' else
                   [(1, column,
                     64 if shortLength else 24 + (column % 2) * 6,
                     63 + shortLength if shortLength else
                     187 - (column % 2) * 6)
                   for column in range(2, 10)]
                   if path != '--synthetic-triangle' else [
            (1, 3, 97, 138), (1, 2, 109, 126),
            (1, 5, 73, 162), (1, 4, 85, 150),
            (1, 7, 48, 187), (1, 6, 60, 175),
            (1, 9, 24, 211), (1, 8, 36, 199),
            (2, 2, 109, 126), (2, 3, 97, 138),
            (2, 4, 85, 150), (2, 5, 73, 162),
            (2, 6, 60, 175), (2, 7, 48, 187),
            (2, 8, 36, 199), (2, 9, 24, 211),
        ])
        for row, col, low, high in records:
            part.createVirtualHelix(row, col, useUndoStack=False)
            part.virtualHelixAtCoord((row, col)).scaffoldStrandSet().createStrand(
                low, high, useUndoStack=False)
        if not os.environ.get('AUTOCS_TEST_REGULAR'):
            part.setGuidedDuplexRegions(dict(
                ((row, col), [(low, high)])
                for row, col, low, high in records))
        document.setSelectedPart(part)
    else:
        with open(path, encoding='utf-8') as stream:
            source = stream.read()
        # Headless legacy decoding otherwise defaults to honeycomb. Respect
        # the explicit extension field used by the modified application.
        sourceObject = json.loads(source)
        helixRange = os.environ.get('AUTOCS_TEST_HELIX_RANGE')
        if helixRange:
            firstNumber, lastNumber = [
                int(value) for value in helixRange.split('-', 1)]
            sourceObject['vstrands'] = [
                record for record in sourceObject['vstrands']
                if firstNumber <= int(record['num']) <= lastNumber]
            for record in sourceObject['vstrands']:
                number = int(record['num'])
                length = len(record['scaf'])
                scaffold = []
                for index in range(length):
                    if number % 2 == 0:
                        scaffold.append([
                            number if index else -1,
                            index - 1 if index else -1,
                            number if index + 1 < length else -1,
                            index + 1 if index + 1 < length else -1])
                    else:
                        scaffold.append([
                            number if index + 1 < length else -1,
                            index + 1 if index + 1 < length else -1,
                            number if index else -1,
                            index - 1 if index else -1])
                record['scaf'] = scaffold
                record['stap'] = [[-1, -1, -1, -1]
                                  for unused in range(length)]
                record['stap_colors'] = []
            source = json.dumps(sourceObject)
        addBases = int(os.environ.get('AUTOCS_TEST_ADD_BASES', '0'))
        if addBases:
            newLength = int(sourceObject.get(
                'num_bases', len(sourceObject['vstrands'][0]['scaf']))) + \
                addBases
            sourceObject['num_bases'] = newLength
            for record in sourceObject['vstrands']:
                number = int(record['num'])
                scaffold = []
                for index in range(newLength):
                    if number % 2 == 0:
                        scaffold.append([
                            number if index else -1,
                            index - 1 if index else -1,
                            number if index + 1 < newLength else -1,
                            index + 1 if index + 1 < newLength else -1])
                    else:
                        scaffold.append([
                            number if index + 1 < newLength else -1,
                            index + 1 if index + 1 < newLength else -1,
                            number if index else -1,
                            index - 1 if index else -1])
                record['scaf'] = scaffold
                record['stap'] = [[-1, -1, -1, -1]
                                  for unused in range(newLength)]
                record['loop'] = [0] * newLength
                record['skip'] = [0] * newLength
                record['scafLoop'] = []
                record['stapLoop'] = []
                record['stap_colors'] = []
            source = json.dumps(sourceObject)
        if sourceObject.get('lattice') == 'square':
            from cadnano2.model.enum import LatticeType
            from cadnano2.model.io.legacydecoder import import_legacy_dict
            import_legacy_dict(document, sourceObject, LatticeType.Square,
                               forceLatticeType=True)
        else:
            decode(document, source)
    part = document.selectedPart()
    if os.environ.get('AUTOCS_TEST_OLD_PROBE'):
        existing = list(part._existingScaffoldCrossoverRecords())
        for fromNumber, toNumber, idx5p, unused_idx3p in existing:
            fromHelix = part.virtualHelix(fromNumber)
            strand5p = fromHelix.scaffoldStrandSet().getStrand(idx5p)
            strand3p = strand5p.connection3p()
            if strand3p is not None and \
                    strand3p.virtualHelix().number() == toNumber:
                part.removeXover(strand5p, strand3p)
        created = part.autoScaffoldCrossovers()
        scaffold = [oligo for oligo in part.oligos()
                    if not oligo.isStaple()]
        print(json.dumps({
            'created': created,
            'components': len(scaffold),
            'loops': sum(1 for oligo in scaffold if oligo.isLoop()),
            'lengths': sorted(oligo.length() for oligo in scaffold),
            'xovers': part._existingScaffoldCrossoverRecords(),
        }), flush=True)
        return
    if os.environ.get('AUTOCS_TEST_SHIFT'):
        shift = int(os.environ['AUTOCS_TEST_SHIFT'])
        for vh in part.getVirtualHelices():
            for strand in list(vh.scaffoldStrandSet()):
                low, high = strand.idxs()
                strand.resize((low + shift, high + shift),
                              useUndoStack=False)
    from cadnano2.model.parts.part import _endpointCrossoverCandidates
    from cadnano2.model.parts.part import _autoScaffoldSnakePaths
    from cadnano2.model.parts.part import _autoScaffoldStraightPaths
    from cadnano2.model.parts.part import _connectAutoScaffoldSnakePaths
    from cadnano2.model.enum import StrandType
    diagnosticHelices = []
    for vh in part.getVirtualHelices():
        neighbors = [neighbor.number() for neighbor in
                     part.getVirtualHelixNeighbors(vh)
                     if neighbor is not None]
        diagnosticHelices.append(
            (vh.number(), vh.coord()[0], vh.coord()[1], neighbors))
    diagnosticPaths = _autoScaffoldSnakePaths(diagnosticHelices)
    diagnosticConnectedSnakePaths = _connectAutoScaffoldSnakePaths(
        diagnosticPaths, diagnosticHelices)
    diagnosticStraightPaths = _autoScaffoldStraightPaths(diagnosticHelices)
    diagnosticConnectedStraightPaths = _connectAutoScaffoldSnakePaths(
        diagnosticStraightPaths, diagnosticHelices)
    diagnosticPotential = sorted(set(
        (vh.number(), neighbor.number(), index)
        for vh in part.getVirtualHelices()
        for neighbor, index, strandType, unused_isLow in
        part.potentialCrossoverList(vh)
        if strandType == StrandType.Scaffold))
    initialStrands = sorted(
        (vh.number(), [strand.idxs() for strand in vh.scaffoldStrandSet()])
        for vh in part.getVirtualHelices())
    initialExistingXovers = part._existingScaffoldCrossoverRecords()
    initialScaffoldOligos = []
    for oligo in part.oligos():
        if oligo.isStaple():
            continue
        strands = list(oligo.strand5p().generator3pStrand())
        strand3p = strands[-1]
        initialScaffoldOligos.append(
            (oligo.strand5p().virtualHelix().number(),
             oligo.strand5p().idx5Prime(),
             strand3p.virtualHelix().number(),
             strand3p.idx3Prime(),
             oligo.length(), oligo.isLoop()))
    initialScaffoldOligos.sort()
    initialOligos, initialCandidates = _endpointCrossoverCandidates(part)
    candidateSummary = sorted(
        (oligo.strand5p().virtualHelix().number(),
         oligo.strand5p().idx5Prime(),
         [(record[2].strand5p().virtualHelix().number(), record[3],
           record[0]) for record in records])
        for oligo, records in initialCandidates.items())
    densityMultiple = int(os.environ.get('AUTOCS_DENSITY_MULTIPLE', '1'))
    if os.environ.get('AUTOCS_TEST_PRE_GUIDED'):
        from cadnano2.model.parts.part import \
            _applyScaffoldCrossoverRecordSet
        _applyScaffoldCrossoverRecordSet(part, set())
        created = part._autoScaffoldCrossoversPreGuided(
            densitySpacing=part._step,
            explicitPaths=diagnosticPaths)
        result = {
            'success': bool(created), 'created': created,
            'message': 'pre-guided probe'}
    elif os.environ.get('AUTOCS_TEST_PREBUILD_ROUTE'):
        part.autoScaffoldCrossovers(
            minimumIndex=None, densityMultiple=densityMultiple,
            rebuildExisting=True, routeOnly=True, returnDetails=True,
            _suppressPreferenceSearch=True)
    if os.environ.get('AUTOCS_TEST_ROUTE_THEN_DENSITY'):
        part.autoScaffoldCrossovers(
            minimumIndex=(int(os.environ['AUTOCS_MINIMUM_INDEX'])
                          if 'AUTOCS_MINIMUM_INDEX' in os.environ else None),
            densityMultiple=densityMultiple, rebuildExisting=True,
            routeOnly=True, returnDetails=True,
            _suppressPreferenceSearch=True)
    if os.environ.get('AUTOCS_TEST_PRE_GUIDED'):
        pass
    elif os.environ.get('AUTOCS_ANALYZE_ONLY'):
        result = {'success': True, 'created': 0, 'message': 'analyze only'}
    else:
        result = part.autoScaffoldCrossovers(
            minimumIndex=(int(os.environ['AUTOCS_MINIMUM_INDEX'])
                          if 'AUTOCS_MINIMUM_INDEX' in os.environ else None),
            densityMultiple=densityMultiple,
            rebuildExisting=(
                not bool(os.environ.get(
                    'AUTOCS_TEST_ROUTE_THEN_DENSITY')) and
                not bool(os.environ.get(
                    'AUTOCS_PRESERVE_EXISTING'))),
            routeOnly=bool(os.environ.get('AUTOCS_ROUTE_ONLY')),
            returnDetails=True,
            selectionVariant=int(os.environ.get('AUTOCS_VARIANT', '0')),
            _suppressPreferenceSearch=bool(
                os.environ.get('AUTOCS_TEST_SUPPRESS_PREFERENCE')),
            _phasePreferenceOffset=(
                int(os.environ['AUTOCS_PHASE_OFFSET'])
                if 'AUTOCS_PHASE_OFFSET' in os.environ else None),
            _unifiedInternal=bool(
                os.environ.get('AUTOCS_TEST_UNIFIED_INTERNAL')))
        if not result.get('success') and \
                os.environ.get('AUTOCS_TEST_REDO_FAILED') and \
                part.undoStack().canRedo():
            part.undoStack().redo()
        if result.get('success') and \
                os.environ.get('AUTOCS_TEST_THIN_SPACING'):
            from cadnano2.model.parts.part import \
                _removeAutoCrossoversOutsideScaffoldLoops, \
                _thinClosedScaffoldLoopToSpacing
            _removeAutoCrossoversOutsideScaffoldLoops(part)
            _thinClosedScaffoldLoopToSpacing(
                part, int(os.environ['AUTOCS_TEST_THIN_SPACING']),
                part._step)
        if result.get('success') and \
                os.environ.get('AUTOCS_TEST_DENSIFY_TARGET'):
            from cadnano2.model.enum import StrandType
            from cadnano2.model.parts.part import \
                _densifyClosedScaffoldLoop
            groupByHelix = dict(
                (number, groupIndex)
                for groupIndex, path in enumerate(diagnosticPaths)
                for number in path)
            routePairs = set(
                tuple(sorted(pair))
                for path in diagnosticPaths
                for pair in zip(path, path[1:]))
            bridgePairs = set(
                tuple(sorted(record[:2]))
                for record in part._existingScaffoldCrossoverRecords()
                if groupByHelix.get(record[0]) !=
                   groupByHelix.get(record[1]))
            finalSeam = set([
                tuple(sorted(diagnosticPaths[-1][-2:]))])
            densePairs = routePairs.difference(finalSeam).union(bridgePairs)
            densityCandidates = []
            seen = set()
            for vh in part.getVirtualHelices():
                strandSet = vh.scaffoldStrandSet()
                is5to3 = strandSet.isDrawn5to3()
                for neighbor, index, strandType, isLowIdx in \
                        part.potentialCrossoverList(vh):
                    if strandType != StrandType.Scaffold:
                        continue
                    fromHelixIs5p = ((isLowIdx and is5to3) or
                                     (not isLowIdx and not is5to3))
                    record = (vh.number(), neighbor.number(),
                              index, index)
                    if not fromHelixIs5p or record in seen or \
                            tuple(sorted(record[:2])) not in densePairs:
                        continue
                    seen.add(record)
                    densityCandidates.append(record)
            _densifyClosedScaffoldLoop(
                part, densityCandidates, part._step,
                6 if part._step == 21 else 7)
            if os.environ.get('AUTOCS_TEST_DENSIFY_DEBUG'):
                from cadnano2.model.parts.part import \
                    _selectAutoScaffoldCrossoverRecords
                print('DENSIFY_DEBUG', json.dumps({
                    'candidate_count': len(densityCandidates),
                    'selected': _selectAutoScaffoldCrossoverRecords(
                        densityCandidates,
                        part._existingScaffoldCrossoverRecords(),
                        part._step,
                        6 if part._step == 21 else 7),
                }), flush=True)
        if result.get('success') and \
                os.environ.get('AUTOCS_TEST_PRUNE_ROWS'):
            from cadnano2.model.parts.part import \
                _pruneAlternatingSquareRowBridges
            _pruneAlternatingSquareRowBridges(
                part, diagnosticStraightPaths)
        removePairs = os.environ.get('AUTOCS_TEST_REMOVE_PAIRS')
        if removePairs:
            from cadnano2.model.parts.part import \
                _mergeRemovedScaffoldXoverBoundaries
            blocked = set(tuple(sorted(int(item) for item in token.split('-')))
                          for token in removePairs.split(','))
            removedBoundaries = set()
            removableEndpoints = set(
                (number, index)
                for fromNumber, toNumber, idx5p, idx3p in
                part._existingScaffoldCrossoverRecords()
                if tuple(sorted((fromNumber, toNumber))) in blocked
                for number, index in (
                    (fromNumber, idx5p), (toNumber, idx3p)))
            for vh in part.getVirtualHelices():
                strands = sorted(list(vh.scaffoldStrandSet()),
                                 key=lambda strand: strand.lowIdx())
                for left, right in zip(strands, strands[1:]):
                    if left.highIdx() + 1 == right.lowIdx() and \
                            ((vh.number(), left.highIdx()) in
                             removableEndpoints or
                             (vh.number(), right.lowIdx()) in
                             removableEndpoints):
                        removedBoundaries.add(
                            (vh.number(), left.highIdx(), right.lowIdx()))
            for fromNumber, toNumber, idx5p, unused_idx3p in list(
                    part._existingScaffoldCrossoverRecords()):
                if tuple(sorted((fromNumber, toNumber))) not in blocked:
                    continue
                fromHelix = part.virtualHelix(fromNumber)
                strand5p = fromHelix.scaffoldStrandSet().getStrand(idx5p)
                strand3p = strand5p.connection3p()
                if strand3p is not None and \
                        strand3p.virtualHelix().number() == toNumber:
                    part.removeXover(strand5p, strand3p)
            _mergeRemovedScaffoldXoverBoundaries(
                part, removedBoundaries)
        if os.environ.get('AUTOCS_TEST_REMERGE_ALLOWED'):
            from cadnano2.model.parts.part import \
                _mergeClosedScaffoldLoops
            allowedPairs = set(
                tuple(sorted(pair))
                for path in diagnosticPaths
                for pair in zip(path, path[1:]))
            _mergeClosedScaffoldLoops(
                part, part._step, 6 if part._step == 21 else 7,
                allowedPairs=allowedPairs)
            if os.environ.get('AUTOCS_TEST_REMERGE_GROUPS'):
                _mergeClosedScaffoldLoops(
                    part, part._step, 6 if part._step == 21 else 7)
    scaffold = [oligo for oligo in part.oligos() if not oligo.isStaple()]
    bridgeTargetApplied = None
    bridgeTarget = os.environ.get('AUTOCS_TEST_BRIDGE_TARGET')
    if bridgeTarget:
        from cadnano2.model.parts.part import \
            _applyScaffoldCrossoverRecordSet
        targetPair = tuple(sorted(
            int(value) for value in
            os.environ['AUTOCS_TEST_BRIDGE_PAIR'].split('-')))
        targetRecords = set(
            part._existingScaffoldCrossoverRecords())
        targetRecords = set(
            record for record in targetRecords
            if tuple(sorted(record[:2])) != targetPair)
        for token in bridgeTarget.split(','):
            fromNumber, toNumber, index = [
                int(value) for value in token.split('-')]
            targetRecords.add(
                (fromNumber, toNumber, index, index))
        bridgeTargetApplied = _applyScaffoldCrossoverRecordSet(
            part, targetRecords)
        scaffold = [
            oligo for oligo in part.oligos()
            if not oligo.isStaple()]
    from cadnano2.model.enum import LatticeType
    from cadnano2.model.io.legacyencoder import legacy_dict_from_part
    from cadnano2.model.io.legacydecoder import import_legacy_dict
    encoded = legacy_dict_from_part(part, 'roundtrip.json')
    roundtripDocument = Document()
    roundtripPart = import_legacy_dict(
        roundtripDocument, encoded,
        LatticeType.Honeycomb if part._step == 21 else LatticeType.Square,
        forceLatticeType=True)
    scaffoldCount = len(scaffold)
    loopCount = sum(1 for oligo in scaffold if oligo.isLoop())
    scaffoldLengths = sorted(oligo.length() for oligo in scaffold)
    scaffoldXovers = part._existingScaffoldCrossoverRecords()
    loopXovers = []
    openXovers = []
    for record in scaffoldXovers:
        vh = part.virtualHelix(record[0])
        strand = vh.scaffoldStrandSet().getStrand(record[2])
        (loopXovers if strand.oligo().isLoop() else openXovers).append(record)
    if os.environ.get('AUTOCS_ASSERT_ADJACENT_ONLY'):
        eligible = set(
            tuple(pair) for pair in result.get('eligible_pairs', []))
        assert all(
            tuple(sorted(record[:2])) in eligible
            for record in scaffoldXovers)
        spacing = int(result['spacing'])
        byDirection = {}
        for record in scaffoldXovers:
            byDirection.setdefault(record[:2], []).append(record[2])
        assert all(
            right - left >= spacing
            for indices in byDirection.values()
            for left, right in zip(sorted(indices), sorted(indices)[1:]))
        directionSpacing = 6 if part._step == 21 else 7
        events = {}
        for fromNumber, toNumber, idx5p, idx3p in scaffoldXovers:
            events.setdefault(fromNumber, []).append((idx5p, toNumber))
            events.setdefault(toNumber, []).append((idx3p, fromNumber))
        for helixEvents in events.values():
            for index, neighbor in helixEvents:
                assert all(
                    otherIndex == index or
                    otherNeighbor == neighbor or
                    abs(otherIndex - index) >= directionSpacing
                    for otherIndex, otherNeighbor in helixEvents)
        if part._step == 32 and \
                os.environ.get('AUTOCS_ASSERT_NO_X8'):
            assert all(record[2] % 8 for record in scaffoldXovers)
        print('ADJACENT_ONLY_ASSERTIONS_OK', flush=True)
    loopRoutes = []
    if os.environ.get('AUTOCS_TEST_LOOP_ROUTES'):
        for oligo in scaffold:
            if not oligo.isLoop():
                continue
            loopRoutes.append([
                (strand.virtualHelix().number(),
                 strand.lowIdx(), strand.highIdx())
                for strand in oligo.strand5p().generator3pStrand()])
    scaffoldOnlyRecords = part.scaffoldOnlyRegionRecords()
    finalPotentialScaffoldXovers = sorted(set(
        (vh.number(), neighbor.number(), index)
        for vh in part.getVirtualHelices()
        for neighbor, index, strandType, unused_isLow in
        part.potentialCrossoverList(vh)
        if strandType == StrandType.Scaffold))
    displayScaffoldOnly = dict(
        (str(vh.coord()), part.scaffoldOnlyIntervals(vh))
        for vh in part.getVirtualHelices())
    stapleTargets = dict(
        (str(vh.coord()), part.stapleTargetSegments(vh))
        for vh in part.getVirtualHelices())
    undoState = None
    if result.get('success') and os.environ.get('AUTOCS_TEST_UNDO'):
        part.undoStack().undo()
        undoState = {
            'scaffold_only_regions': part.scaffoldOnlyRegionRecords(),
            'scaffold_xovers':
                part._existingScaffoldCrossoverRecords(),
            'scaffold_oligos': len([
                oligo for oligo in part.oligos()
                if not oligo.isStaple()]),
            'loops': sum(1 for oligo in part.oligos()
                         if not oligo.isStaple() and oligo.isLoop()),
        }
        if os.environ.get('AUTOCS_ASSERT_UNDO'):
            assert undoState['scaffold_xovers'] == initialExistingXovers
            print('AUTOCS_UNDO_ASSERTIONS_OK', flush=True)
    if os.environ.get('AUTOCS_COMPACT'):
        print(json.dumps({
            'variant': int(os.environ.get('AUTOCS_VARIANT', '0')),
            'result': result,
            'loops': loopCount,
            'lengths': scaffoldLengths,
            'scaffold_xovers': scaffoldXovers,
            'loop_xovers': loopXovers,
            'open_xovers': openXovers,
            'bridge_target_applied': bridgeTargetApplied,
            'final_potential_scaffold_xovers':
                finalPotentialScaffoldXovers,
            'loop_routes': loopRoutes,
            'diagnostic_paths': diagnosticPaths,
            'diagnostic_helices': diagnosticHelices,
            'diagnostic_connected_snake_paths':
                diagnosticConnectedSnakePaths,
            'diagnostic_straight_paths': diagnosticStraightPaths,
            'diagnostic_connected_straight_paths':
                diagnosticConnectedStraightPaths,
            'diagnostic_directions': dict(
                (vh.number(),
                 vh.scaffoldStrandSet().isDrawn5to3())
                for vh in part.getVirtualHelices()),
            'undo_state': undoState,
        }, ensure_ascii=False), flush=True)
        return
    print(json.dumps({
        'result': result,
        'scaffold_oligos': scaffoldCount,
        'loops': loopCount,
        'lengths': scaffoldLengths,
        'scaffold_xovers': scaffoldXovers,
        'scaffold_only_regions': scaffoldOnlyRecords,
        'display_scaffold_only': displayScaffoldOnly,
        'roundtrip_scaffold_only_regions':
            roundtripPart.scaffoldOnlyRegionRecords(),
        'undo_state': undoState,
        'staple_target_segments': stapleTargets,
        'initial_endpoint_oligos': len(initialOligos),
        'initial_endpoint_candidates': candidateSummary,
        'initial_scaffold_strands': initialStrands,
        'initial_scaffold_xovers': initialExistingXovers,
        'initial_scaffold_oligos': initialScaffoldOligos,
        'panel_bounds': [part.minBaseIdx(), part.maxBaseIdx()],
        'diagnostic_paths': diagnosticPaths,
        'diagnostic_potential_scaffold_xovers': diagnosticPotential,
    }, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
